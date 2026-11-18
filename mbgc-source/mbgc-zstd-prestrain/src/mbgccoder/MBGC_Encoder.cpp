#include "MBGC_Encoder.h"

#include <numeric>

#include "../matching/input_with_libdeflate_wrapper.h"
#include "../matching/SimpleSequenceMatcher.h"

#include "../coders/PropsLibrary.h"

#include <fstream>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <omp.h>

MBGC_Encoder::MBGC_Encoder(MBGC_Params *mbgcParams):
    MultipleGenomeMatchingProcessor(mbgcParams), params(mbgcParams) {

}

void MBGC_Encoder::print_invalid_kseq_status_message(int kseq_status)
{
    fprintf(stderr, "Consider using lossy compression mode (-%c option).\n", MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT);
}

void MBGC_Encoder::setProteinsProfile()
{
    params->setProteinsCompressionProfile();
}

void MBGC_Encoder::initStreamsForG0Ref() {
    targetLiterals.clear();
    fileHeaders.clear();
    fileHeadersTemplates.clear();
    targetLiterals.resize(1);
    fileHeaders.push_back("");
    fileHeadersTemplates.push_back("");
}

void MBGC_Encoder::processG0RefContig(const char *seq, size_t len) {
    targetLiterals[0].append(seq, len);
    targetLiterals[0].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
}

void MBGC_Encoder::updateHeadersTemplate(uint32_t fileIndex, const string &header) {
    string& currentHeadersTemplate = fileHeadersTemplates[fileIndex];
    if (currentHeadersTemplate.empty()) {
        currentHeadersTemplate = header;
        return;
    }
    int hPos = 0;
    int templatePos = currentHeadersTemplate[0] == MBGC_Params::MATCH_MARK ? 1 : 0;
    string newTemplate(currentHeadersTemplate, 0, templatePos);
    while (templatePos + MBGC_Params::MINIMAL_PATTERN_LENGTH <= currentHeadersTemplate.size()) {
        int patternEnd = currentHeadersTemplate.find(MBGC_Params::MATCH_MARK, templatePos);
        if (patternEnd == std::string::npos)
            patternEnd = currentHeadersTemplate.size();
        if (patternEnd - templatePos < MBGC_Params::MINIMAL_PATTERN_LENGTH)
            templatePos = patternEnd + 1;
        else {
            string part = currentHeadersTemplate.substr(templatePos + MBGC_Params::MAX_PATTERN_SHIFT,
                                                        MBGC_Params::MINIMAL_PATTERN_LENGTH - MBGC_Params::MAX_PATTERN_SHIFT);
            int tmpH = header.find(part, hPos);
            if (tmpH != std::string::npos) {
                int tmpT = templatePos + MBGC_Params::MAX_PATTERN_SHIFT;
                while (tmpH != hPos && tmpT != templatePos
                       && currentHeadersTemplate[--tmpT] == header[--tmpH]);
                templatePos = tmpT + (currentHeadersTemplate[tmpT] == header[tmpH] ? 0 : 1);
                if (tmpH != hPos && newTemplate.empty())
                    newTemplate.push_back(MBGC_Params::MATCH_MARK);
                while (tmpH != header.size() && tmpT != patternEnd
                       && currentHeadersTemplate[++tmpT] == header[++tmpH]);
                if (tmpT - templatePos >= MBGC_Params::MINIMAL_PATTERN_LENGTH) {
                    newTemplate.append(currentHeadersTemplate.substr(templatePos, tmpT - templatePos));
                    if (tmpH != header.size() || tmpT != currentHeadersTemplate.size())
                        newTemplate.push_back(MBGC_Params::MATCH_MARK);
                }
                templatePos = tmpT + 1;
                hPos = tmpH + 1;
            } else {
                templatePos += MBGC_Params::MAX_PATTERN_SHIFT;
            }
        }
    }
    if (newTemplate.back() != MBGC_Params::MATCH_MARK && hPos < header.size())
        newTemplate.push_back(MBGC_Params::MATCH_MARK);
    currentHeadersTemplate = newTemplate;
}

void MBGC_Encoder::processHeader(uint32_t fileIndex, const string& header) {
    updateHeadersTemplate(fileIndex, header);
    fileHeaders[fileIndex].append(header);
    fileHeaders[fileIndex].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
}

void MBGC_Encoder::processFileName(string &fileName) {
    switch (fileName[0]) {
    case ' ': case '\t': case '\r': case '\n':
        fprintf(stderr, "ERROR: Sequences list file %s is invalid!\n", params->seqListFileName.c_str());
        fprintf(stderr, "File name cannot start from '%c' whitespace: %s\n", fileName[0],
            fileName.c_str());
        exit(EXIT_FAILURE);
    default: ;
    }
    for (char c : fileName) {
        switch (c) {
        case '<': case '>': case ':': case '"': case '|': case '?': case '*':
            fprintf(stderr, "ERROR: Sequences list file %s is invalid!\n", params->seqListFileName.c_str());
            fprintf(stderr, "It contains invalid file name: %s\n", fileName.c_str());
            exit(EXIT_FAILURE);
        default: ;
        }
    }
    int pos, len;
    PgHelpers::normalizePath(fileName, pos, len, params->ignoreFastaFilesPath);
    if (fileNamesSet.emplace(fileName.substr(pos, len)).second == true) {
        namesStr.append(fileName, pos, len);
        namesStr.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
    } else {
        fprintf(stderr, "WARNING: skipped %s file (already in the archive).\n", fileName.substr(pos, len).c_str());
        fileName.clear();
    }
}

void MBGC_Encoder::processTargetMeta(uint32_t seqCount, uint64_t dnaLineLength) {
    PgHelpers::writeValue<uint32_t>(seqsCountsDest, seqCount);
    PgHelpers::writeValue<uint64_t>(dnaLineLengthsDest, dnaLineLength);
}

size_t extensionsMatchedCharsAll = 0;
size_t extensionsMismatchesAll = 0;
size_t removedGapBreakingMatchesAll = 0;

void MBGC_Encoder::processLiteral(char *destPtr, uint32_t pos, uint64_t length, size_t destLen, string &refExtRes) {
    if (params->isLiteralProperForRefExtension(length)) {
        int before = params->refLiteralBeforeAfterExt < pos ? params->refLiteralBeforeAfterExt : pos;
        int after = (pos + length + params->refLiteralBeforeAfterExt < destLen) ?
                    params->refLiteralBeforeAfterExt : destLen - pos - length;
        refExtRes.append(destPtr + pos - before, length + before + after);
    }
}

size_t MBGC_Encoder::getMatchLoadedPos(size_t pos) {
    while (pos + (matcher->getRefLength() - matcher->REF_SHIFT) < refExtLoadedPosArr.back())
        pos += (matcher->getRefLength() - matcher->REF_SHIFT);
    return pos;
}

size_t MBGC_Encoder::processMatches(vector<PgTools::TextMatch>& textMatches, char *destStart, size_t destLen, int targetIdx,
                                    size_t matchingLockPos) {
    uint32_t pos = 0;
    int64_t unmatchedChars = 0;
    uint32_t extensionsMatchedChars = 0;
    uint32_t extensionsMismatches = 0;
    uint32_t totalDestOverlap = 0;
    uint32_t totalMatched = 0;
    uint32_t removedGapBreakingMatches = 0;
    int64_t jj = 0, j;
    for (j = 0; j < textMatches.size(); j++, jj++) {
        TextMatch &match = textMatches[j];
#ifdef DEVELOPER_BUILD
        if (match.posDestText < pos) {
            fprintf(stderr, "WARNING: Unexpected matches state: overlapping matches\n");
            uint32_t overflow = pos - match.posDestText;
            if (overflow >= match.length) {
                totalDestOverlap += match.length;
                jj--;
                continue;
            }
            totalDestOverlap += overflow;
            match.shiftStartPos(overflow);
        }
        if (match.length < params->k) {
            fprintf(stderr, "WARNING: Unexpected matches state: too short match %ld < %d \n", match.length, params->k);
            totalDestOverlap += match.length;
            jj--;
            continue;
        }
#endif
        if (params->enableExtensionsWithMismatches) {
            if (j + 1 < textMatches.size()) {
                if (jj > 0 && textMatches[j + 1].pairedWith(textMatches[jj - 1]) &&
                              !match.pairedWith(textMatches[jj - 1]) &&
                              match.length < params->gapBreakingMatchMinLength) {
                    int64_t leftExtension = 0;
                    if (match.endPosDestText() == textMatches[j + 1].posDestText) {
                        const char *srcPtr = matcher->getRef() + textMatches[j + 1].posSrcText;
                        const char *ptr = destStart + textMatches[j + 1].posDestText;
                        while (*--ptr == *--srcPtr)
                            leftExtension++;
                        textMatches[j + 1].shiftStartPos(-leftExtension);
                    }
                    removedGapBreakingMatches++;
                    jj--;
                    continue;
                }
            }
        }
        totalMatched += match.length;
        uint64_t length = match.posDestText - pos;
        unmatchedChars += length;
        pos = match.endPosDestText();
        textMatches[jj] = textMatches[j];
    }
    textMatches.resize(jj);
    uint64_t literalsLeft = destLen - pos;
    unmatchedChars += literalsLeft;
    int unmatchedFractionFactor = unmatchedFractionFactors[2 * (targetsCountShift + targetIdx)];
    if (processedTargetsCount < targetIdx - params->allowedTargetsOutrunForDissimilarContigs &&
            params->isContigDissimilar(destLen, unmatchedChars, unmatchedFractionFactor))
        return PROCESSING_MATCHES_SKIPPED_DUE_TO_CONTIG_DISSIMILARITY;
    bool pairedGap[MBGC_Params::MAX_GAP_DEPTH] = {};
    const int NO_GAP = -1;
    int64_t gapStartIdx = NO_GAP;
    int64_t gapEndIdx = NO_GAP;
    int gapCurIdx = 0;
    bool loadLiteralsToRef = !params->isContigProperForRefExtension(destLen, unmatchedChars, unmatchedFractionFactor);
    pos = 0;
    bool isGap = false;
    for (int64_t j = 0; j < textMatches.size(); j++) {
        TextMatch &match = textMatches[j];
        literalsLeft = match.posDestText - pos;
        if (params->enableExtensionsWithMismatches) {
            if (!isGap && literalsLeft) {
                literalsLeft -= extendMatchLeft(destStart, literalsLeft, match, targetIdx, matchingLockPos, extensionsMatchedChars,
                                                extensionsMismatches);
            }
            if (j == gapEndIdx) {
                gapStartIdx = NO_GAP;
                gapEndIdx = NO_GAP;
            }
        }
        targetLiterals[targetIdx].append(destStart + pos, literalsLeft);
        if (loadLiteralsToRef) processLiteral(destStart, pos, literalsLeft, destLen, targetRefExtensions[targetIdx]);
        bool skipOffset = pairedGap[gapCurIdx];
        if (!skipOffset) {
            PgHelpers::writeValue<uint32_t>(targetMapOffDests[targetIdx], match.posSrcText);
            if (params->enable40bitReference)
                targetMapOff5thByte[targetIdx].push_back(match.posSrcText >> 32);
        }
        if (params->frugal64bitLenEncoding)
            PgHelpers::writeUInt64Frugal(targetMapLenDests[targetIdx], match.length);
        else
            PgHelpers::writeValue<uint32_t>(targetMapLenDests[targetIdx], match.length);
        targetLiterals[targetIdx].push_back(MBGC_Params::MATCH_MARK);
        pos = match.endPosDestText();

        literalsLeft = (j + 1 < textMatches.size() ? textMatches[j + 1].posDestText : destLen) - pos;
        int gapIdx = (gapCurIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;
        int gCnt = std::min((int) (textMatches.size() - j - 1), (int) params->gapDepthOffsetEncoding);
        int reduce = 0, g;
        for (g = 1; g <= gCnt; g++, gapIdx = (gapIdx + 1) % MBGC_Params::MAX_GAP_DEPTH) {
            if (pairedGap[gapIdx] || (!params->lazyDecompressionSupport && g == 1 && gapStartIdx == NO_GAP && literalsLeft == 0)) {
                reduce++;
                continue;
            }
            if (match.pairedWith(textMatches[j + g], matchingLockPos)) {
                if (params->lazyDecompressionSupport) {
                    if(!match.nextSrcRegionLoadingPos) {
                        auto itr = std::upper_bound(refExtLoadedPosArr.begin(), refExtLoadedPosArr.end(),
                                                    getMatchLoadedPos(match.posSrcText));
                        match.nextSrcRegionLoadingPos = itr == refExtLoadedPosArr.end() ? UINT64_MAX : *itr;
                    }
                    if (getMatchLoadedPos(textMatches[j + g].posSrcText) >= match.nextSrcRegionLoadingPos)
                        continue;
                    textMatches[j + g].nextSrcRegionLoadingPos = match.nextSrcRegionLoadingPos;
                }
                pairedGap[gapIdx] = true;
                if (params->enableExtensionsWithMismatches && gapEndIdx <= j + g && g <= params->gapDepthMismatchesEncoding) {
                    gapStartIdx = j;
                    gapEndIdx = j + g;
                }
                break;
            }
        }
        if (gCnt)
            targetGapDeltas[targetIdx].push_back(g <= gCnt ? g - reduce : 0);
        pairedGap[gapCurIdx] = false;
        gapCurIdx = (gapCurIdx + 1) % MBGC_Params::MAX_GAP_DEPTH;

        bool gapStart = gapStartIdx == j;
        bool gapEnd = gapEndIdx == j + 1;
        bool gapMiddle = gapStartIdx < j && j + 1 < gapEndIdx;
        isGap = gapStart || gapMiddle || gapEnd;
        if (params->enableExtensionsWithMismatches) {
            TextMatch& corrMatch = isGap ? textMatches[gapStartIdx] : match;
            uint64_t extensionLength = extendMatchRight(destStart + pos, corrMatch, match,
                                                        literalsLeft, targetIdx, isGap, gapStart, gapMiddle, gapEnd,
                                                        extensionsMatchedChars, extensionsMismatches);
            literalsLeft -= extensionLength;
            pos += extensionLength;
        }
    }
    literalsLeft = destLen - pos;
    targetLiterals[targetIdx].append(destStart + pos, literalsLeft);
    if (loadLiteralsToRef) processLiteral(destStart, pos, literalsLeft, destLen, targetRefExtensions[targetIdx]);
    textMatches.clear();

#pragma omp atomic update
    unmatchedCharsAll += unmatchedChars;
#pragma omp atomic update
    extensionsMatchedCharsAll += extensionsMatchedChars;
#pragma omp atomic update
    extensionsMismatchesAll += extensionsMismatches;
#pragma omp atomic update
    totalMatchedAll += totalMatched;
#pragma omp atomic update
    totalDestOverlapAll += totalDestOverlap;
#pragma omp atomic update
    totalDestLenAll += destLen;
#pragma omp atomic update
    removedGapBreakingMatchesAll += removedGapBreakingMatches;
    return unmatchedChars;
}

uint64_t MBGC_Encoder::extendMatchRight(const char *gapStartPtr, const TextMatch &coreMatch, const TextMatch &match, uint64_t &length,
                                        int targetIdx, bool isGap, bool gapStart, bool gapMiddle, bool gapEnd,
                                        uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches) {
    if (length == 0) {
        if (gapMiddle)
            targetGapMismatchesFlags[targetIdx].push_back(1);
        return 0;
    }
    const char *srcStartPtr = matcher->getRef();
    const char *gapPtr = gapStartPtr;
    const char *srcPtr = matcher->getRef() + (isGap
            ? coreMatch.posSrcText + match.endPosDestText() - coreMatch.posDestText
            : match.endPosSrcText());
    size_t tmpPos = matcher->getLoadingPosition();
    const char *srcLoadingPos = matcher->getRef() + tmpPos;
    const char *srcGuard = srcPtr + length;
    const char *validSrcGuard = srcPtr + length;
    if (srcPtr == srcLoadingPos)
        validSrcGuard = srcPtr;
    if (!isGap) {
        const char *srcEnd = matcher->getRef() + matcher->getMaxRefLength();
        if (validSrcGuard > srcEnd)
            validSrcGuard = srcEnd;
        if (srcPtr <= srcLoadingPos && srcLoadingPos < validSrcGuard)
            validSrcGuard = srcLoadingPos;
    }
    if (gapStart || !isGap) {
        if (params->lazyDecompressionSupport && *srcPtr == MBGC_Params::REF_REGION_SEPARATOR)
            validSrcGuard = srcPtr;
        extensionsMismatches++;
        targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion && srcPtr < validSrcGuard ?
                                            mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
        gapPtr++;
    } else
        srcPtr--;
    int mismatches_score = params->mmsMismatchesInitialScore;
    while (++srcPtr < validSrcGuard && (!params->lazyDecompressionSupport || *srcPtr != MBGC_Params::REF_REGION_SEPARATOR) &&
        (isGap || mismatches_score < params->mmsMismatchesScoreThreshold)) {
        bool mismatch = *gapPtr != *srcPtr;
        targetGapMismatchesFlags[targetIdx].push_back(mismatch ? 1 : 0);
        if (mismatch) {
            extensionsMismatches++;
            mismatches_score += params->mmsMismatchPenalty;
            targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion ?
                    mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
        } else {
            extensionsMatchedChars++;
            mismatches_score -= params->mmsMatchBonus;
            if (mismatches_score < 0) mismatches_score = 0;
        }
        gapPtr++;
    }
    while (srcPtr++ < srcGuard && (isGap || mismatches_score < params->mmsMismatchesScoreThreshold)) {
        targetGapMismatchesFlags[targetIdx].push_back(1);
        targetLiterals[targetIdx].push_back(*gapPtr++);
        extensionsMismatches++;
        mismatches_score += params->mmsMismatchPenalty;
    }
    if ((isGap && !gapEnd) || (!isGap && mismatches_score < params->mmsMismatchesScoreThreshold))
        targetGapMismatchesFlags[targetIdx].push_back(1);
    return gapPtr - gapStartPtr;
}

uint64_t MBGC_Encoder::extendMatchLeft(const char *destStart, uint64_t length, const TextMatch &match,
                                   int targetIdx, size_t matchingLockPos,
                                   uint32_t &extensionsMatchedChars, uint32_t &extensionsMismatches) {
    const char *srcMatch = matcher->getRef() + match.posSrcText;
    const char *srcGuard = matcher->getRef() + 1;
    if (srcGuard < srcMatch - MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH)
        srcGuard = srcMatch - MBGC_Params::MAX_EXTEND_MATCH_LEFT_LENGTH;
    const char *srcLock = matcher->getRef() + matchingLockPos;
    if (srcGuard < srcLock && srcLock <= srcMatch)
        srcGuard = srcLock;
    bool isGuardKnownInDecompression = true;
    if (srcGuard < srcMatch - length) {
        srcGuard = srcMatch - length;
        isGuardKnownInDecompression = false;
    }
    if (srcGuard == srcMatch)
        return 0;
    const char *srcPtr = srcMatch - 1;
    const char* gapPtr = destStart + match.posDestText - 1;
    bool validSrcRegion = !params->lazyDecompressionSupport || *srcPtr != MBGC_Params::REF_REGION_SEPARATOR;
    extensionsMismatches++;
    targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion && validSrcRegion ?
            mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
    int mismatches_score = params->mmsMismatchesInitialScore;
    while (validSrcRegion && srcPtr > srcGuard && mismatches_score < params->mmsMismatchesScoreThreshold) {
        bool mismatch = *--gapPtr != *--srcPtr;
        if (params->lazyDecompressionSupport && *srcPtr == MBGC_Params::REF_REGION_SEPARATOR) {
            validSrcRegion = false;
            srcPtr++;
            gapPtr++;
            break;
        }
        targetGapMismatchesFlags[targetIdx].push_back(mismatch ? 1 : 0);
        if (mismatch) {
            mismatches_score += params->mmsMismatchPenalty;
            extensionsMismatches++;
            targetLiterals[targetIdx].push_back(params->mismatchesWithExclusion ?
                    mismatchesCoder->mismatch2code(*srcPtr, *gapPtr) : *gapPtr);
        } else {
            mismatches_score -= params->mmsMatchBonus;
            if (mismatches_score < 0) mismatches_score = 0;
            extensionsMatchedChars++;
        }
    }
    while (!validSrcRegion && srcPtr > srcGuard && mismatches_score < params->mmsMismatchesScoreThreshold) {
        srcPtr--;
        targetGapMismatchesFlags[targetIdx].push_back(1);
        targetLiterals[targetIdx].push_back(*--gapPtr);
        extensionsMismatches++;
        mismatches_score += params->mmsMismatchPenalty;
    }
    if ((srcPtr != srcGuard || !isGuardKnownInDecompression) && mismatches_score < params->mmsMismatchesScoreThreshold)
        targetGapMismatchesFlags[targetIdx].push_back(1);
    return srcMatch - srcPtr;
}

void MBGC_Encoder::loadFileNames() {
    if (!params->repackCommand) {
        if (singleFastaFileMode) {
            fileNames.push_back(params->inputFileName);
            if (params->appendCommand) {
                filesCount = 1;
                return;
            }
        } else {
            ifstream listSrc(params->seqListFileName, ios_base::in | ios_base::binary);
            if (listSrc.fail()) {
                fprintf(stderr, "cannot open sequences list file %s\n", params->seqListFileName.c_str());
                exit(EXIT_FAILURE);
            }
            string line;
            while (getline(listSrc, line))
                fileNames.push_back(line);
            if (fileNames.empty()) {
                fprintf(stderr, "ERROR: Sequences list file %s is empty!\n", params->seqListFileName.c_str());
                exit(EXIT_FAILURE);
            }
            if (fileNames[0][0] == '>') {
                fprintf(stderr, "ERROR: Sequences list file %s is invalid!\n", params->seqListFileName.c_str());
                fprintf(stderr, "File name cannot start from '>' character: %s\n", fileNames[0].c_str());
                fprintf(stderr, "TIP: To compress a single FASTA file use -%c option.\n",
                    MBGC_Params::INPUT_FILE_OPT);
                exit(EXIT_FAILURE);
            }
        }
    }

    MultipleGenomeMatchingProcessor::loadFileNames();

    if (appendedFilesCount + filesCount == 1) {
        fprintf(stderr, "Switching to single file compression mode.\n");
        singleFastaFileMode = true;
        params->lazyDecompressionSupport = false;
    }
}

void MBGC_Encoder::initProcessTargetsWithParallelIO() {
    fileHeaders.resize(filesCount);
    fileHeadersTemplates.resize(filesCount);
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        vector<ostringstream> &dests = p.second;
        dests.resize(1);
    }
    for(vector<string>* s: usedByteStreams)
        s->resize(1);
}

void MBGC_Encoder::finalizeProcessTargetsWithParallelIO() {
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        string &stream = p.first;
        vector<ostringstream> &dests = p.second;
        stream = dests[0].str();
        dests.pop_back();
    }
}

void MBGC_Encoder::processAfterSequence(uint32_t targetIdx) {
    targetLiterals[targetIdx].push_back(MBGC_Params::SEQ_SEPARATOR_MARK);
}

void MBGC_Encoder::processAfterTarget(uint32_t targetIdx) {
    if (params->enableExtensionsWithMismatches)
        targetGapMismatchesFlags[targetIdx].push_back(MBGC_Params::FILE_SEPARATOR_MARK);
}

void MBGC_Encoder::processAfterTargetWithParallelIO(size_t matcherLoaderStartPos) {
    processAfterTarget(0);
    if (params->lazyDecompressionSupport) {
        matcher->loadSeparator(MBGC_Params::REF_REGION_SEPARATOR);
        size_t refExtSize = matcher->getLoadedRefLength() - matcherLoaderStartPos;
        PgHelpers::writeUInt64Frugal(refExtSizeDest, refExtSize);
        refExtLoadedPosArr.emplace_back(refExtLoadedPosArr.back() + refExtSize);
    }
    size_t tmp = matcher->acquireWorkerMatchingLockPos();
    locksPosStream.append((char*) &tmp, sizeof(tmp));
    matcher->releaseWorkerMatchingLockPos(tmp);
}

void MBGC_Encoder::initProcessTarget(uint32_t targetIdx) {
    targetLiterals[targetIdx].reserve(refG0InitPos / 32);
}


void MBGC_Encoder::initParallelProcessing() {
    MultipleGenomeMatchingProcessor::initParallelProcessing();
    fileHeaders.resize(singleFastaFileMode ? targetsCount + targetsShift : filesCount);
    fileHeadersTemplates.resize(singleFastaFileMode ? targetsCount + targetsShift : filesCount);
    for(vector<string>* s: usedByteStreams)
        s->resize(targetsCount);
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        vector<ostringstream> &dests = p.second;
        dests.resize(targetsCount);
    }
    if (params->lazyDecompressionSupport)
        refExtLoadedPosArr.reserve(targetsCountShift + targetsCount + targetsShift);
}

void MBGC_Encoder::finalizeParallelProcessingInSingleFastaFileMode() {
    MultipleGenomeMatchingProcessor::finalizeParallelProcessingInSingleFastaFileMode();
    fileHeaders.resize(targetsCount + targetsShift);
    fileHeadersTemplates.resize(targetsCount + targetsShift);
    for(pair<string&, vector<ostringstream>&> p: usedStreamsWithDests) {
        vector<ostringstream> &dests = p.second;
        dests.resize(targetsCount);
    }
    for(vector<string>* s: usedByteStreams)
        s->resize(targetsCount);
}

void MBGC_Encoder::finalizeParallelProcessingOfTarget(uint32_t targetIdx, size_t matcherLoaderStartPos) {
    for (vector<string> *s: usedByteStreams) {
        if (targetIdx) {
            (*s)[0].append((*s)[targetIdx]);
            (*s)[targetIdx].clear();
            (*s)[targetIdx].shrink_to_fit();
        }
    }
    for (pair<string &, vector<ostringstream> &> p: usedStreamsWithDests) {
        string &stream = p.first;
        vector<ostringstream> &dests = p.second;
        stream.append(dests[targetIdx].str());
        dests[targetIdx].str("");
        dests[targetIdx].clear();
    }
    if (params->lazyDecompressionSupport) {
        matcher->loadSeparator(MBGC_Params::REF_REGION_SEPARATOR);
        size_t refExtSize = matcher->getLoadedRefLength() - matcherLoaderStartPos;
        PgHelpers::writeUInt64Frugal(refExtSizeDest, refExtSize);
        refExtLoadedPosArr.emplace_back(matcher->getLoadedRefLength());
    }
    locksPosStream.append((char *) &matchingLocksPos[targetIdx], sizeof(matchingLocksPos[targetIdx]));
};

void MBGC_Encoder::applyTemplateToHeaders(uint32_t fileIdx) {
    vector<size_t> matchesPos;
    string& header = fileHeaders[fileIdx];
    string& hTemplate = fileHeadersTemplates[fileIdx];
    matchesPos.clear();
    size_t matchPos = 0;
    while ((matchPos = hTemplate.find(MBGC_Params::MATCH_MARK, matchPos)) != std::string::npos)
        matchesPos.push_back(matchPos++);
    size_t headersPos = 0;
    size_t newPos = 0;
    matchesPos.push_back(hTemplate.size());
    while (headersPos < header.size()) {
        int tPos = 0;
        size_t headerEnd = header.find(MBGC_Params::SEQ_SEPARATOR_MARK, headersPos);
        for (int i = 0; i < matchesPos.size() - 1; i++) {
            headersPos += matchesPos[i] - tPos;
            tPos = matchesPos[i] + 1;
            int hPos = matchesPos[i + 1] - tPos == 0 ? headerEnd :
                       header.find(hTemplate.substr(tPos, matchesPos[i + 1] - tPos), headersPos);
            std::copy(header.begin() + headersPos, header.begin() + hPos, header.begin() + newPos);
            newPos += hPos - headersPos;
            headersPos = hPos;
            header[newPos++] = MBGC_Params::MATCH_MARK;
        }
        headersPos = headerEnd + 1;
    }
    header.resize(newPos);
}

void MBGC_Encoder::prepareHeadersStreams() {
    int64_t count = singleFastaFileMode && !params->sequentialMatching ? targetsCount + targetsShift : filesCount;
    auto stringSize = [](int a, const std::string& b) { return a + b.size(); };
    size_t totalHeadersLength = std::accumulate(std::next(fileHeaders.begin()), fileHeaders.end(),
                                                fileHeaders[0].size(), stringSize);
    size_t totalHeaderTemplatesLength = std::accumulate(std::next(fileHeadersTemplates.begin()),
                                                        fileHeadersTemplates.end(), fileHeadersTemplates[0].size(), stringSize);
    headersStr.reserve(totalHeadersLength + count);
    headersTemplates.reserve(totalHeaderTemplatesLength + count);
    for(uint32_t f = 0; f < count; f++) {
        headersStr.append(fileHeaders[f]);
        headersTemplates.append(fileHeadersTemplates[f]);
        if (params->storeFileSeparatorMarksInHeadersStream)
            headersStr.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
        headersTemplates.push_back(MBGC_Params::FILE_SEPARATOR_MARK);
    }
    fileHeaders.clear();
    fileHeadersTemplates.clear();
}

void MBGC_Encoder::prepareAndCompressStreams() {
    string seqsCount = seqsCountsDest.str();
    seqsCountsDest.clear();
    string dnaLineLengthStr = dnaLineLengthsDest.str();
    dnaLineLengthsDest.clear();
    ostream* out;
    string tempArchiveFileName = params->outArchiveFileName + MBGC_Params::TEMPORARY_FILE_SUFFIX;
    if (params->isStdoutMode(false)) {
        out = &cout;
    } else {
        PgHelpers::createFolders(tempArchiveFileName);
        out = new fstream(tempArchiveFileName, ios::out | ios::binary | ios::trunc);
    }
    if (!*out)
    {
        fprintf(stderr, "ERROR: Cannot create archive: %s!\n", params->outArchiveFileName.c_str());
        exit(EXIT_FAILURE);
    }
    params->write(*out);
    writeStats(*out);
    prepareHeadersStreams();
    string rcMapOff, rcMapLen;
    if (params->rcRedundancyRemoval)
        SimpleSequenceMatcher::rcMatchSequence(targetLiterals[0], rcMapOff, rcMapLen, params->rcMatchMinLength);

    params->enableOmpThreads(params->backendThreads);
    vector<CompressionJob> cJobs;
    bool fastDecoder = params->coderMode == MBGC_Params::SPEED_CODER_MODE;
    bool repoMode = params->coderMode == MBGC_Params::REPO_CODER_MODE;
    bool bestRatio = params->coderMode >= MBGC_Params::MAX_CODER_MODE;
    auto namesCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 6);
    cJobs.emplace_back("file names stream... ", namesStr, namesCoderProps.get());
    auto seqCountersCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("sequences counters stream... ", seqsCount, seqCountersCoderProps.get());
    auto headersTemplateCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    cJobs.emplace_back("sequences headers templates stream... ", headersTemplates, headersTemplateCoderProps.get());
    auto headers1stCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_NORMAL, LZMA_DATAPERIODCODE_8_t);
    auto headers2ndCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_NORMAL, 3);
    auto headersCoderProps = params->ultraStreamsCompression ?
            getCompoundCoderProps(headers1stCoderProps.get(), headers2ndCoderProps.get()) :
            getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,  fastDecoder ? 8 : 16);
    int headersBlocksCount = (bestRatio ? 1 : 2) ;
    ParallelBlocksCoderProps blockHeadersCoderProps(headersBlocksCount, headersCoderProps.get());
    cJobs.emplace_back("sequences headers stream... ", headersStr, &blockHeadersCoderProps);
    auto dnaLineLengthsProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 8);
    cJobs.emplace_back("dna line lengths stream... ", dnaLineLengthStr, dnaLineLengthsProps.get());
    auto refExtFactorsProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    cJobs.emplace_back("unmatched fraction factors stream... ", unmatchedFractionFactors.data(),
                                   unmatchedFractionFactors.size(), refExtFactorsProps.get());
    auto seqCoderProps = params->ultraStreamsCompression ?
            getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_MAX, LZMA_DATAPERIODCODE_8_t)
            :(fastDecoder || params->k == MBGC_Params::PROTEINS_PROFILE_KMER_LENGTH ? getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_FAST, LZMA_DATAPERIODCODE_8_t)
                : getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, fastDecoder ? 4 :
                                       (params->enableExtensionsWithMismatches ?
                                           (params->mismatchesWithExclusion ? 5 : 7) :
                                    (params->sequentialMatching ? 6 : 7))));
    int seqBlocksCount = bestRatio ? 1 : (fastDecoder ? 4 : (repoMode ? 3 : 2));
    ParallelBlocksCoderProps blockSeqCoderProps(seqBlocksCount, seqCoderProps.get());
    cJobs.emplace_back("literals stream... ", targetLiterals[0], &blockSeqCoderProps);
    auto rcMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_MAX, LZMA_DATAPERIODCODE_32_t);
    auto rcMapLenCoderProps = getDefaultCoderProps(LZMA_CODER, CODER_LEVEL_MAX, LZMA_DATAPERIODCODE_8_t);
    if (params->rcRedundancyRemoval) {
        cJobs.emplace_back("rc matches offsets stream... ", rcMapOff, rcMapOffCoderProps.get());
        cJobs.emplace_back("rc matches lengths stream... ", rcMapLen, rcMapLenCoderProps.get());
    }
    auto refLocksCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 8);
    cJobs.emplace_back("reference loading locks position stream... ", locksPosStream, refLocksCoderProps.get());
    auto gapDeltaCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                                   bestRatio ? 12 : (fastDecoder ? 2 : 8));
    int gapDeltaBlocksCount = repoMode ? 4 : (bestRatio ? 1 : 2);
    ParallelBlocksCoderProps blockGapDeltaCoderProps(gapDeltaBlocksCount, gapDeltaCoderProps.get());
    cJobs.emplace_back("gap between paired matches (fwd-delta) stream... ", targetGapDeltas[0], &blockGapDeltaCoderProps);
    auto gapMismatchesFlagsCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX,
                                                             fastDecoder ? 8 : 14);
    int gapMismatchesFlagsBlocksCount = fastDecoder ? 3 : 2;
    ParallelBlocksCoderProps blockGapMismatchesFlagsCoderProps(gapMismatchesFlagsBlocksCount, gapMismatchesFlagsCoderProps.get());
    uint8_t defaultCoderLevel = bestRatio ? CODER_LEVEL_MAX : (fastDecoder ? CODER_LEVEL_FAST : CODER_LEVEL_NORMAL);
    if (params->enableExtensionsWithMismatches)
        cJobs.emplace_back("mismatches in gaps between matches (flags) stream... ", targetGapMismatchesFlags[0], &blockGapMismatchesFlagsCoderProps);
    auto nMapOffCoderProps = getDefaultCoderProps(LZMA_CODER, defaultCoderLevel, LZMA_DATAPERIODCODE_32_t);
    int mapOffBlocksCount = (bestRatio ? 1 : (fastDecoder ? 4 : 3));
    ParallelBlocksCoderProps blockNMapOffCoderProps(mapOffBlocksCount, nMapOffCoderProps.get());
    cJobs.emplace_back("matches offsets stream... ", mapOffStream, &blockNMapOffCoderProps);
    auto mapOff5thByteCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 4);
    if (refFinalTotalLength > UINT32_MAX)
        cJobs.emplace_back("matches offsets 5th byte stream... ", targetMapOff5thByte[0], mapOff5thByteCoderProps.get());
    auto nMapLenCoderProps = getDefaultCoderProps(LZMA_CODER, defaultCoderLevel,
        params->frugal64bitLenEncoding ? LZMA_DATAPERIODCODE_16_t : LZMA_DATAPERIODCODE_32_t);
    int mapLenBlocksCount = (bestRatio ? 2 : (fastDecoder ? 10 : (repoMode ? 3 : 5)));
    ParallelBlocksCoderProps blockNMapLenCoderProps(mapLenBlocksCount, nMapLenCoderProps.get());
    cJobs.emplace_back("matches lengths stream... ", mapLenStream, &blockNMapLenCoderProps);
    string refExtSizeStream = refExtSizeDest.str();
    auto refExtSizeCoderProps = getDefaultCoderProps(PPMD7_CODER, CODER_LEVEL_MAX, 6);
    if (params->lazyDecompressionSupport)
        cJobs.emplace_back("reference extensions size stream... ", refExtSizeStream, refExtSizeCoderProps.get());
    CompressionJob::writeCompressedCollectiveParallel(*out, cJobs);
    size_t outSize = out->tellp();
    if (!params->isStdoutMode(false)) {
        delete (out);
        if (std::ifstream(params->outArchiveFileName)) {
            if (!params->appendCommand && !params->forceOverwrite) {
                fprintf(stderr, "ERROR: file %s already exists (use -%c to force overwrite)\n",
                    params->outArchiveFileName.data(), MBGC_Params::FORCE_OVERWRITE_OPT);
                exit(EXIT_FAILURE);
            }
            remove(params->outArchiveFileName.c_str());
        }
        if (rename(tempArchiveFileName.c_str(), params->outArchiveFileName.c_str()) != 0) {
            fprintf(stderr, "Error preparing output file: %s\n", params->outArchiveFileName.c_str());
            exit(EXIT_FAILURE);
        };
    }
    *PgHelpers::devout << endl;
    *PgHelpers::appout << "compressed " << totalFilesLength << " bytes to " << PgHelpers::getPartShareString(outSize, totalFilesLength) << endl;
    *PgHelpers::logout << PgHelpers::time_millis() << "       \t";
    *PgHelpers::logout << PgHelpers::toString((totalFilesLength * 1.0 / outSize + 0.05), 1) << "\t";
    *PgHelpers::logout << totalFilesLength << "\t" << outSize << "\t";
}

void MBGC_Encoder::writeStats(ostream &out) const {
    PgHelpers::writeValue(out, appendedFilesCount + filesCount);
    PgHelpers::writeValue(out, totalFilesLength);
    PgHelpers::writeValue(out, refG0InitPos);
    PgHelpers::writeValue(out, largestFileLength);
    PgHelpers::writeValue(out, largestContigSize);
    PgHelpers::writeValue(out, refFinalTotalLength);
    PgHelpers::writeValue(out, refBuffersCount);
    PgHelpers::writeValue(out, (uint64_t) targetLiterals[0].size());
}

void MBGC_Encoder::printAdditionalMatchingStats() {
    if (params->enableExtensionsWithMismatches) {
        *PgHelpers::devout << "removed matches breaking gaps total: " << removedGapBreakingMatchesAll << endl;
        *PgHelpers::devout << "extensions matched chars: " << extensionsMatchedCharsAll << endl;
        *PgHelpers::devout << "extensions mismatched chars: ~" <<
                           PgHelpers::getPartShareString(extensionsMismatchesAll,
                                             extensionsMatchedCharsAll + extensionsMismatchesAll) << endl;
        *PgHelpers::devout << "final unmatched chars: " << (unmatchedCharsAll - extensionsMatchedCharsAll)
                           << " (removed: " <<
                           PgHelpers::getPartShareString(totalMatchedAll + extensionsMatchedCharsAll, totalDestLenAll) << ")"
                           << endl;
    }
}

void MBGC_Encoder::encode() {
    if (!params->noOfThreadsLimited) {
        params->coderThreads = omp_get_num_procs();
        params->backendThreads = omp_get_num_procs();
    }
    params->enableOmpThreads(params->coderThreads);
    if (!params->isStdoutMode(false) && !params->forceOverwrite
        && std::ifstream(params->outArchiveFileName)
        && (!params->appendCommand || params->inArchiveFileName != params->outArchiveFileName)) {
        fprintf(stderr, "ERROR: file %s already exists (use -%c to force overwrite)\n",
                    params->outArchiveFileName.data(), MBGC_Params::FORCE_OVERWRITE_OPT);
        exit(EXIT_FAILURE);
    }
    if (PgHelpers::numberOfThreads == 1)
        params->setSequentialMatchingMode();
    if (!params->appendCommand) {
        singleFastaFileMode = !params->inputFileName.empty();
        loadFileNames();
        loadG0Ref(fileNames[0]);
    }
    enrollByteStream(targetLiterals);
    enrollStream(mapOffStream, targetMapOffDests);
    if (params->enable40bitReference)
        enrollByteStream(targetMapOff5thByte);
    enrollStream(mapLenStream, targetMapLenDests);
    enrollByteStream(targetGapDeltas);
    if (params->enableExtensionsWithMismatches) {
        enrollByteStream(targetGapMismatchesFlags);
        params->initMismatchesMatchingScoreParams();
    }
    if (params->lazyDecompressionSupport) {
        refExtLoadedPosArr.emplace_back(matcher->getLoadingPosition());
    }
    if (params->currentUnmatchedFractionFactor < params->unmatchedFractionRCFactor) {
        fprintf(stderr, "WARNING: %c cannot be smaller than %c. Decreased %c down to %d.\n",
                MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT, MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT,
                MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT, params->currentUnmatchedFractionFactor);
        params->unmatchedFractionRCFactor = params->currentUnmatchedFractionFactor;
    }
    if (params->matcherWorkingThreads > PgHelpers::numberOfThreads) {
        fprintf(stderr, "WARNING: %c cannot be larger than number of threads limit. Decreased %c down to %d.\n",
                MBGC_Params::MATCHER_NO_OF_THREADS_OPT,
                MBGC_Params::MATCHER_NO_OF_THREADS_OPT, PgHelpers::numberOfThreads);
        params->matcherWorkingThreads = PgHelpers::numberOfThreads;
    }
    performMatching();

    prepareAndCompressStreams();
}

template<bool lazyMode>
void MBGC_Encoder::appendTemplate(MBGC_Decoder<lazyMode>* decoder, MBGC_Params& newParams) {
    decoder->decode();
    namesStr = std::move(decoder->namesStr);
    size_t namesPos = 0;
    fileNamesSet.reserve(decoder->filesCount * 2);
    for (int i = 0; i < decoder->filesCount; i++) {
        size_t tmpNamesPos = namesStr.find(MBGC_Params::FILE_SEPARATOR_MARK, namesPos);
        string name = namesStr.substr(namesPos, tmpNamesPos - namesPos);
        fileNamesSet.emplace(name);
        namesPos = tmpNamesPos + 1;
    }
    seqsCountsDest << decoder->seqsCount;
    headersTemplates = std::move(decoder->headersTemplates);
    headersStr = std::move(decoder->headersStr);
    dnaLineLengthsDest << decoder->dnaLineLengthsStr;
    unmatchedFractionFactors.assign(decoder->unmatchedFractionFactorsArray,
                                    decoder->unmatchedFractionFactorsArray + decoder->unmatchedFractionFactorsStr.size());
    params->currentUnmatchedFractionFactor = *((uint8_t*) unmatchedFractionFactors.data() + unmatchedFractionFactors.size() - 2);
    params->unmatchedFractionRCFactor = *((uint8_t*) unmatchedFractionFactors.data() + unmatchedFractionFactors.size() - 1);
    targetLiterals.push_back(std::move(decoder->literalStr));
    locksPosStream = std::move(decoder->matchingLocksPosStream);
    targetGapDeltas.push_back(std::move(decoder->gapDeltaStr));
    if (params->enableExtensionsWithMismatches)
        targetGapMismatchesFlags.push_back(std::move(decoder->gapMismatchesFlagsStr));
    targetMapOffDests.resize(1);
    targetMapOffDests[0] << decoder->mapOffStream;
    if (decoder->refTotalLength > UINT32_MAX)
        targetMapOff5thByte.push_back(std::move(decoder->mapOff5thByteStream));
    targetMapLenDests.resize(1);
    targetMapLenDests[0] << decoder->mapLenStream;
    if (params->lazyDecompressionSupport) {
        refExtSizeDest << decoder->refExtSizeStream;
        refExtLoadedPosArr.assign(decoder->refExtLoadedPosArr.begin(),
                                  decoder->refExtLoadedPosArr.end() - 1);
    }

    singleFastaFileMode = !params->inputFileName.empty();
    appendedFilesCount = singleFastaFileMode ? 0 : decoder->filesCount;
    targetsCountShift = decoder->targetsCount;
    totalFilesLength = decoder->totalFilesLength;
    refG0InitPos = decoder->refG0InitPos;
    largestFileLength = decoder->largestFileLength;
    largestContigSize = decoder->largestContigLength;
    refBuffersCount = decoder->refBuffersCount;

    loadFileNames();
    string refStr = std::move(decoder->refStr);
    size_t initialRefLength = refStr.size();
    size_t fileSize = openFirstFile(fileNames[0].c_str(), true);
    if (decoder->reachedRefLengthCount == 0) {
        int64_t targetsProportion = singleFastaFileMode ?
                                (((int64_t) fileSize) - 1) / ((int64_t) initialRefLength) + 1 :
                                filesCount / appendedFilesCount + 1;
        int tmp = 15 - (__builtin_clz(targetsProportion) / 3);
        tmp = tmp < 5 ? 5 : (tmp > 12 ? 12 : tmp);
        params->referenceFactor = (size_t) 1 << tmp;
    } else
        params->referenceFactor = 1;
    if (newParams.coderThreads > 0)
        params->setThreadsLimit(newParams.coderThreads);
    if (newParams.matcherWorkingThreads > 0)
        params->setMatcherThreads(newParams.matcherWorkingThreads);
    if (newParams.sequentialMatching || params->g0IsTarget || params->matcherWorkingThreads == 1)
        params->setSequentialMatchingMode();
    if (!params->sequentialMatching && singleFastaFileMode)
        switchToSplitIterMode();
    size_t basicRefLength = refG0InitPos;
    basicRefLength = basicRefLength > params->MIN_BASIC_BLOCK_SIZE ? basicRefLength : params->MIN_BASIC_BLOCK_SIZE;
    targetsCount = singleFastaFileMode && !params->sequentialMatching ?
        (((int64_t) fileSize) + basicRefLength - 1) / basicRefLength : filesCount;
    targetsShift = 0;
    initMatcher(refStr.data() + 1, refStr.size() - 1, initialRefLength, true);
    if (decoder->reachedRefLengthCount > 0) {
        matcher->setPosition(decoder->refPos, decoder->reachedRefLengthCount);
    }
    if (basicRefLength > params->AVG_REF_INIT_SIZE * 2)
        readingBufferSize = 1 + readingBufferSize * params->AVG_REF_INIT_SIZE * 2 / basicRefLength;
    *PgHelpers::appout << "processed reference - " << PgHelpers::time_millis() << " [ms]" << endl;
    delete(decoder);
    encode();
}

void MBGC_Encoder::append(MBGC_Params& newParams) {
    MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(params);
    params->engageAppendCommand();
    if (params->isLazyDecompressionEnabled())
        appendTemplate<true>((MBGC_Decoder<true>*) decoder, newParams);
    else
        appendTemplate<false>((MBGC_Decoder<false>*) decoder, newParams);
}

template<bool lazyMode>
void MBGC_Encoder::repackTemplate(MBGC_Decoder<lazyMode>* decoder) {
    decoder->params->repackCommand = true;
    iterableDecoder = decoder;
    decoder->iterateInit();
    decoder->getSelectedFilesNames(fileNames);
    encode();
}

void MBGC_Encoder::repack(MBGC_Params& inParams) {
    MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&inParams);
    if (inParams.isLazyDecompressionEnabled())
        repackTemplate<true>((MBGC_Decoder<true>*) decoder);
    else
        repackTemplate<false>((MBGC_Decoder<false>*) decoder);
}


// ************************ G0 site-wise comparison (streaming, no decompression to disk) ************************

// Each FASTA record is treated as one strain; the first record of the first file is G0.
static void parseFastaRecords(const uint8_t* data, size_t size, vector<pair<string, string>>& records) {
    size_t i = 0;
    while (i < size) {
        while (i < size && data[i] != '>')
            i++;
        if (i >= size)
            break;
        i++;
        const size_t hStart = i;
        while (i < size && data[i] != '\n' && data[i] != '\r')
            i++;
        string name((const char*) data + hStart, i - hStart);
        const size_t sp = name.find_first_of(" \t"); // keep the first token as strain id
        if (sp != string::npos)
            name = name.substr(0, sp);
        while (i < size && (data[i] == '\n' || data[i] == '\r'))
            i++;
        string seq;
        while (i < size && data[i] != '>') {
            char c = (char) data[i++];
            if (c == '\n' || c == '\r' || c == ' ' || c == '\t')
                continue;
            if (c >= 'a' && c <= 'z')
                c -= 32; // uppercase
            seq.push_back(c);
        }
        records.emplace_back(move(name), move(seq));
    }
}

static inline bool isComparableBase(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

// LEB128 varint helpers for the compact diff-site dump
static void writeVarint(ostream& out, uint64_t v) {
    while (v >= 0x80) {
        out.put((char) (v | 0x80));
        v >>= 7;
    }
    out.put((char) v);
}

static bool readVarint(const uint8_t*& p, const uint8_t* end, uint64_t& v) {
    v = 0;
    int shift = 0;
    while (p < end) {
        uint8_t b = *p++;
        v |= (uint64_t) (b & 0x7F) << shift;
        if (!(b & 0x80))
            return true;
        shift += 7;
        if (shift > 63)
            return false;
    }
    return false;
}

template<typename T>
static void writeLE(ostream& out, T v) {
    for (size_t i = 0; i < sizeof(T); i++) {
        out.put((char) (v & 0xFF));
        v >>= 8;
    }
}

template<typename T>
static bool readLE(const uint8_t*& p, const uint8_t* end, T& v) {
    if ((size_t) (end - p) < sizeof(T))
        return false;
    v = 0;
    for (size_t i = 0; i < sizeof(T); i++)
        v |= (T) *p++ << (8 * i);
    return true;
}

static void writePosList(ostream& out, const vector<uint64_t>& pos) {
    writeLE<uint32_t>(out, (uint32_t) pos.size());
    uint64_t prev = 0;
    for (size_t i = 0; i < pos.size(); i++) {
        writeVarint(out, pos[i] - prev); // delta coding (first entry stores the absolute position)
        prev = pos[i];
    }
}

// load 1-based site positions, return sorted unique 0-based positions
static vector<uint64_t> loadSitesFile(const string& sitesFileName) {
    ifstream sitesIn(sitesFileName, ios::in | ios::binary | ios::ate);
    if (!sitesIn) {
        fprintf(stderr, "ERROR: cannot open sites file %s\n", sitesFileName.c_str());
        exit(EXIT_FAILURE);
    }
    const size_t fileSize = sitesIn.tellg();
    sitesIn.seekg(0);
    string buf(fileSize, '\0');
    sitesIn.read(&buf[0], fileSize);
    if (!sitesIn && fileSize > 0) {
        fprintf(stderr, "ERROR: failed to read sites file %s\n", sitesFileName.c_str());
        exit(EXIT_FAILURE);
    }
    vector<uint64_t> sites;
    sites.reserve(fileSize / 4 + 1);
    const char* p = buf.data();
    const char* end = p + buf.size();
    uint64_t lineNo = 1;
    auto isSpace = [](char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; };
    while (p < end) {
        while (p < end && isSpace(*p)) {
            lineNo += *p == '\n';
            p++;
        }
        if (p >= end)
            break;
        if (*p == '#') { // comment line
            while (p < end && *p != '\n')
                p++;
            continue;
        }
        const char* tokStart = p;
        uint64_t v = 0;
        while (p < end && *p >= '0' && *p <= '9')
            v = v * 10 + (uint64_t) (*p++ - '0');
        if (p == tokStart || (p < end && !isSpace(*p)) || v == 0) {
            const char* tokEnd = p;
            while (tokEnd < end && !isSpace(*tokEnd))
                tokEnd++;
            fprintf(stderr, "ERROR: invalid 1-based site position \"%.*s\" in %s (line %llu).\n",
                    (int) (tokEnd - tokStart), tokStart, sitesFileName.c_str(), (unsigned long long) lineNo);
            exit(EXIT_FAILURE);
        }
        sites.push_back(v - 1); // convert to 0-based
    }
    if (sites.empty()) {
        fprintf(stderr, "ERROR: no sites found in %s\n", sitesFileName.c_str());
        exit(EXIT_FAILURE);
    }
    sort(sites.begin(), sites.end());
    sites.erase(unique(sites.begin(), sites.end()), sites.end());
    return sites;
}

template<bool lazyMode>
void MBGC_Encoder::compareTemplate(MBGC_Decoder<lazyMode>* decoder, const vector<uint64_t>& sites, ostream& out,
                                   ostream* dumpOut) {
    decoder->iterateInit();
    vector<string> names;
    decoder->getSelectedFilesNames(names);
    if (names.empty()) {
        fprintf(stderr, "ERROR: no files found in the archive.\n");
        exit(EXIT_FAILURE);
    }
    string g0seq;
    uint64_t strainCount = 0;
    uint64_t differentLengthStrains = 0;
    bool dumpHeaderWritten = false;
    for (size_t i = 0; i < names.size(); i++) {
        uint8_t* res = nullptr;
        size_t size = 0;
        if (!decoder->iterateNext(names[i].c_str(), res, size)) {
            fprintf(stderr, "ERROR: failed to decode file %s from the archive.\n", names[i].c_str());
            exit(EXIT_FAILURE);
        }
        vector<pair<string, string>> records;
        parseFastaRecords(res, size, records);
        free(res);
        if (records.empty())
            fprintf(stderr, "WARNING: no FASTA records found in file %s, skipped.\n", names[i].c_str());
        for (auto& rec : records) {
            if (strainCount == 0) {
                g0seq = rec.second;
                if (!sites.empty() && sites.back() >= g0seq.size()) {
                    fprintf(stderr, "ERROR: site %llu exceeds the first (G0) record sequence length (%llu bp).\n",
                            (unsigned long long) (sites.back() + 1), (unsigned long long) g0seq.size());
                    exit(EXIT_FAILURE);
                }
                out << "#g0_record\t" << rec.first << "\n#sites_total\t" << sites.size()
                    << "\n#strain\tdiff_count\tvalid_sites\n";
            }
            strainCount++;
            const string& seq = rec.second;
            if (seq.size() != g0seq.size())
                differentLengthStrains++;
            const size_t limit = min(seq.size(), g0seq.size());
            const char* a = g0seq.data();
            const char* b = seq.data();
            uint64_t diff = 0, valid = 0;
            vector<uint64_t> diffPos, invalidPos; // used only when dumping
            const bool dump = dumpOut != nullptr;
            for (uint64_t s : sites) {
                if (s >= limit) {
                    if (dump)
                        invalidPos.push_back(s);
                    continue;
                }
                const char ca = a[s], cb = b[s];
                if (!isComparableBase(ca) || !isComparableBase(cb)) {
                    if (dump)
                        invalidPos.push_back(s);
                    continue; // N, gap, etc. - skipped
                }
                valid++;
                if (ca != cb) {
                    diff++;
                    if (dump)
                        diffPos.push_back(s);
                }
            }
            out << rec.first << '\t' << diff << '\t' << valid << '\n';
            if (dump) {
                if (!dumpHeaderWritten) { // written when processing the first (G0) record
                    dumpOut->write("mDIF", 4);
                    writeLE<uint32_t>(*dumpOut, 1); // format version
                    writeLE<uint64_t>(*dumpOut, g0seq.size());
                    writeLE<uint16_t>(*dumpOut, (uint16_t) rec.first.size());
                    dumpOut->write(rec.first.data(), rec.first.size());
                    dumpHeaderWritten = true;
                }
                writeLE<uint16_t>(*dumpOut, (uint16_t) rec.first.size());
                dumpOut->write(rec.first.data(), rec.first.size());
                writePosList(*dumpOut, diffPos);
                writePosList(*dumpOut, invalidPos);
            }
        }
    }
    if (differentLengthStrains > 0)
        fprintf(stderr, "WARNING: %llu record(s) have sequence length different from the first (G0) record; "
                        "out-of-range sites were skipped for them.\n", (unsigned long long) differentLengthStrains);
}

static ostream* openOutputFile(const string& outFileName, bool forceOverwrite, ofstream& outFile) {
    if (outFileName.empty() || outFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS)
        return &cout;
    if (!forceOverwrite && ifstream(outFileName)) {
        fprintf(stderr, "ERROR: file %s already exists (use -%c to force overwrite)\n",
                outFileName.c_str(), MBGC_Params::FORCE_OVERWRITE_OPT);
        exit(EXIT_FAILURE);
    }
    outFile.open(outFileName, ios::out | ios::trunc | ios::binary);
    if (!outFile) {
        fprintf(stderr, "ERROR: cannot create a file %s\n", outFileName.c_str());
        exit(EXIT_FAILURE);
    }
    return &outFile;
}

void MBGC_Encoder::compare(MBGC_Params& inParams, const string& sitesFileName, const string& outFileName,
                           const string& dumpFileName) {
    vector<uint64_t> sites = loadSitesFile(sitesFileName);

    ofstream outFile;
    ostream* out = openOutputFile(outFileName, inParams.forceOverwrite, outFile);

    ofstream dumpFile;
    ostream* dumpOut = nullptr;
    if (!dumpFileName.empty()) {
        dumpOut = openOutputFile(dumpFileName, inParams.forceOverwrite, dumpFile);
    }

    inParams.repackCommand = true; // redirects decoded content to the iterator instead of files
    MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&inParams);
    if (inParams.isLazyDecompressionEnabled())
        compareTemplate<true>((MBGC_Decoder<true>*) decoder, sites, *out, dumpOut);
    else
        compareTemplate<false>((MBGC_Decoder<false>*) decoder, sites, *out, dumpOut);
    delete(decoder);
    if (outFile.is_open())
        outFile.close();
    if (dumpFile.is_open())
        dumpFile.close();
}

// Fast query on a previously dumped diff-site file: no archive access needed.
// The raw diff is loaded into memory exactly once; every (sites, output) job
// then re-parses that resident buffer (no zstd decompression, no disk re-read
// per sample), so a batch of N samples costs one diff load instead of N.

namespace {

struct LoadedDiff {
    string buf;      // whole raw diff (header + record names + varint stream), resident once
    uint64_t g0len = 0;
    string g0name;
    size_t dataOff = 0;  // offset of the first record
};

static bool readWholeFile(const string& fileName, string& buf) {
    ifstream in(fileName, ios::in | ios::binary | ios::ate);
    if (!in)
        return false;
    const size_t fileSize = in.tellg();
    in.seekg(0);
    buf.assign(fileSize, '\0');
    if (fileSize) {
        in.read(&buf[0], fileSize);
        if (!in)
            return false;
    }
    return true;
}

static bool loadDiffFile(const string& diffBinFileName, LoadedDiff& d) {
    if (!readWholeFile(diffBinFileName, d.buf)) {
        fprintf(stderr, "ERROR: cannot open diff file %s\n", diffBinFileName.c_str());
        return false;
    }
    const uint8_t* p = (const uint8_t*) d.buf.data();
    const uint8_t* end = p + d.buf.size();
    if (d.buf.size() < 18 || memcmp(p, "mDIF", 4) != 0) {
        fprintf(stderr, "ERROR: %s is not a valid diff file (bad magic).\n", diffBinFileName.c_str());
        return false;
    }
    p += 4;
    uint32_t version;
    uint64_t g0len;
    uint16_t g0NameLen;
    if (!readLE<uint32_t>(p, end, version) || version != 1 || !readLE<uint64_t>(p, end, g0len)
        || !readLE<uint16_t>(p, end, g0NameLen) || (size_t) (end - p) < g0NameLen) {
        fprintf(stderr, "ERROR: %s is corrupted (bad header).\n", diffBinFileName.c_str());
        return false;
    }
    d.g0len = g0len;
    d.g0name.assign((const char*) p, g0NameLen);
    p += g0NameLen;
    d.dataOff = (size_t) (p - (const uint8_t*) d.buf.data());
    return true;
}

// One (sites, output) job: re-parses the resident raw diff and counts.
// The varint stream cannot be split, so parsing is sequential; the counting
// of each parallel block uses `innerThreads` threads.
static void runDiffQueryJob(const LoadedDiff& d, const string& diffBinFileName,
                            const pair<string, string>& job, bool forceOverwrite, int innerThreads) {
    struct QueryRec {
        uint64_t nameOff;
        uint32_t nameLen;
        uint64_t diffOff, diffCnt, invOff, invCnt;
    };
    const size_t BLOCK_POS_LIMIT = 1u << 24; // ~16M positions (~128 MB) per block

    const string& sitesFileName = job.first;
    const string& outFileName = job.second;
    vector<uint64_t> sites = loadSitesFile(sitesFileName);

    uint64_t outOfRange = 0;
    for (uint64_t s : sites)
        if (s >= d.g0len)
            outOfRange++;
    if (outOfRange > 0)
        fprintf(stderr, "WARNING: %llu site(s) exceed the G0 sequence length (%llu bp) "
                        "and are excluded from valid sites.\n",
                (unsigned long long) outOfRange, (unsigned long long) d.g0len);

    ofstream outFile;
    ostream* out = openOutputFile(outFileName, forceOverwrite, outFile);

    const uint64_t inRangeSites = sites.size() - outOfRange;

    // Site membership bitset over genome positions -> O(1) lookup per diff position
    // (fallback to binary search when the genome is too large for a dense bitset).
    const bool useBitset = d.g0len <= (1ULL << 32);
    vector<uint64_t> siteBits;
    if (useBitset) {
        siteBits.assign((size_t) ((d.g0len + 63) / 64), 0);
        for (uint64_t s : sites)
            if (s < d.g0len)
                siteBits[(size_t) (s >> 6)] |= 1ULL << (s & 63);
    }

    const uint8_t* p = (const uint8_t*) d.buf.data() + d.dataOff;
    const uint8_t* end = (const uint8_t*) d.buf.data() + d.buf.size();

    vector<QueryRec> recs;
    vector<uint64_t> posBuf;
    string outBuf;

    auto flushBlock = [&]() {
        if (recs.empty())
            return;
        vector<string> lines(recs.size());
#pragma omp parallel for num_threads(innerThreads) schedule(dynamic, 32)
        for (int64_t i = 0; i < (int64_t) recs.size(); i++) {
            const QueryRec& r = recs[i];
            const uint64_t* diffPos = posBuf.data() + r.diffOff;
            const uint64_t* invalidPos = posBuf.data() + r.invOff;
            uint64_t diff = 0, invalid = 0;
            if (useBitset) {
                const uint64_t* bits = siteBits.data();
                for (uint64_t k = 0; k < r.diffCnt; k++) {
                    const uint64_t s = diffPos[k];
                    if (s < d.g0len)
                        diff += (bits[(size_t) (s >> 6)] >> (s & 63)) & 1;
                }
                for (uint64_t k = 0; k < r.invCnt; k++) {
                    const uint64_t s = invalidPos[k];
                    if (s < d.g0len)
                        invalid += (bits[(size_t) (s >> 6)] >> (s & 63)) & 1;
                }
            } else {
                for (uint64_t k = 0; k < r.diffCnt; k++)
                    if (diffPos[k] < d.g0len && binary_search(sites.begin(), sites.end(), diffPos[k]))
                        diff++;
                for (uint64_t k = 0; k < r.invCnt; k++)
                    if (invalidPos[k] < d.g0len && binary_search(sites.begin(), sites.end(), invalidPos[k]))
                        invalid++;
            }
            string& line = lines[i];
            line.append(d.buf.data() + r.nameOff, r.nameLen);
            line.push_back('\t');
            line += to_string(diff);
            line.push_back('\t');
            line += to_string(inRangeSites - invalid);
            line.push_back('\n');
        }
        for (const string& line : lines)
            outBuf.append(line);
        recs.clear();
        posBuf.clear();
    };

    outBuf.append("#g0_record\t").append(d.g0name).append("\n#sites_total\t")
          .append(to_string(sites.size())).append("\n#strain\tdiff_count\tvalid_sites\n");

    const uint8_t* bufStart = (const uint8_t*) d.buf.data();
    auto decodePosList = [&](uint64_t& off, uint64_t& cnt) {
        uint32_t n;
        if (!readLE<uint32_t>(p, end, n)) {
            fprintf(stderr, "ERROR: %s is corrupted (truncated position list).\n", diffBinFileName.c_str());
            exit(EXIT_FAILURE);
        }
        off = posBuf.size();
        cnt = n;
        posBuf.resize(posBuf.size() + n);
        uint64_t prev = 0, delta;
        for (uint32_t k = 0; k < n; k++) {
            if (!readVarint(p, end, delta)) {
                fprintf(stderr, "ERROR: %s is corrupted (truncated position list).\n", diffBinFileName.c_str());
                exit(EXIT_FAILURE);
            }
            prev += delta;
            posBuf[off + k] = prev;
        }
    };

    while (p < end) {
        uint16_t nameLen;
        if (!readLE<uint16_t>(p, end, nameLen) || (size_t) (end - p) < nameLen) {
            fprintf(stderr, "ERROR: %s is corrupted (truncated record).\n", diffBinFileName.c_str());
            exit(EXIT_FAILURE);
        }
        QueryRec r;
        r.nameOff = (uint64_t) (p - bufStart);
        r.nameLen = nameLen;
        p += nameLen;
        decodePosList(r.diffOff, r.diffCnt);
        decodePosList(r.invOff, r.invCnt);
        recs.push_back(r);
        if (posBuf.size() >= BLOCK_POS_LIMIT) {
            flushBlock();
            if (outBuf.size() >= (8u << 20)) {
                out->write(outBuf.data(), outBuf.size());
                outBuf.clear();
            }
        }
    }
    flushBlock();
    out->write(outBuf.data(), outBuf.size());
    if (outFile.is_open())
        outFile.close();
}

// Batch (N >= 2) single-pass counting over the union of all samples' sites:
// - one diff load + one varint parse (not N)
// - mask[p] stores which samples contain position p (bit i = sample i)
// - each diff/invalid position updates all owning samples at once
// Output is byte-identical to running each sample separately.
static void runDiffQueryJobsUnion(const LoadedDiff& d, const string& diffBinFileName,
                                  const vector<pair<string, string>>& jobs,
                                  bool forceOverwrite, int noOfThreads) {
    const size_t N = jobs.size();
    if (N < 2)
        return;

    vector<vector<uint64_t>> sitesList;
    sitesList.reserve(N);
    vector<uint64_t> inRangeSites(N, 0);
    for (size_t i = 0; i < N; i++) {
        sitesList.push_back(loadSitesFile(jobs[i].first));
        uint64_t outOfRange = 0;
        for (uint64_t s : sitesList.back())
            if (s >= d.g0len)
                outOfRange++;
        if (outOfRange > 0)
            fprintf(stderr, "WARNING: %llu site(s) in %s exceed the G0 sequence length (%llu bp) "
                            "and are excluded from valid sites.\n",
                    (unsigned long long) outOfRange, jobs[i].first.c_str(),
                    (unsigned long long) d.g0len);
        inRangeSites[i] = sitesList.back().size() - outOfRange;
    }

    // union site mask: bit i of mask[p] = sample i has site p
    // (dense over the genome; one 64-bit word per 64 samples per position)
    const uint32_t maskWords = (uint32_t) ((N + 63) / 64);
    const uint64_t maskEntries = (uint64_t) d.g0len * maskWords;
    if (maskEntries > (1ull << 31)) { // ~16 GB cap for the dense mask
        fprintf(stderr, "ERROR: union mask (%llu entries) too large for genome size %llu bp; "
                        "use single-sample queries instead.\n",
                (unsigned long long) maskEntries, (unsigned long long) d.g0len);
        exit(EXIT_FAILURE);
    }
    vector<uint64_t> mask((size_t) maskEntries, 0);
    for (size_t i = 0; i < N; i++) {
        const uint64_t word = (uint64_t) i >> 6, bit = 1ULL << (i & 63);
        for (uint64_t s : sitesList[i])
            if (s < d.g0len)
                mask[(size_t) s * maskWords + word] |= bit;
    }

    // open all outputs upfront (forceOverwrite semantics preserved)
    vector<ofstream> outFiles(N);
    vector<ostream*> outs(N);
    for (size_t i = 0; i < N; i++)
        outs[i] = openOutputFile(jobs[i].second, forceOverwrite, outFiles[i]);

    vector<string> outBufs(N);
    for (size_t i = 0; i < N; i++) {
        outBufs[i].append("#g0_record\t").append(d.g0name).append("\n#sites_total\t")
                  .append(to_string(sitesList[i].size())).append("\n#strain\tdiff_count\tvalid_sites\n");
    }

    struct QueryRec {
        uint64_t nameOff;
        uint32_t nameLen;
        uint64_t diffOff, diffCnt, invOff, invCnt;
    };
    const size_t BLOCK_POS_LIMIT = 1u << 24; // ~16M positions (~128 MB) per block
    vector<QueryRec> recs;
    vector<uint64_t> posBuf;

    auto flushBlock = [&]() {
        if (recs.empty())
            return;
        // diffMat/invMat flat arrays [rec * N + sample]; each record is owned
        // by exactly one thread -> no races, no atomics.
        vector<uint64_t> diffMat(recs.size() * N, 0), invMat(recs.size() * N, 0);
#pragma omp parallel for num_threads(noOfThreads) schedule(dynamic, 32)
        for (int64_t i = 0; i < (int64_t) recs.size(); i++) {
            const QueryRec& r = recs[i];
            uint64_t* dRow = diffMat.data() + (size_t) i * N;
            uint64_t* iRow = invMat.data() + (size_t) i * N;
            const uint64_t* diffPos = posBuf.data() + r.diffOff;
            const uint64_t* invalidPos = posBuf.data() + r.invOff;
            for (uint64_t k = 0; k < r.diffCnt; k++) {
                const uint64_t s = diffPos[k];
                if (s >= d.g0len)
                    continue;
                const uint64_t* mp = mask.data() + (size_t) s * maskWords;
                for (uint32_t w = 0; w < maskWords; w++) {
                    uint64_t m = mp[w];
                    while (m) {
                        const int b = __builtin_ctzll(m);
                        dRow[(size_t) w * 64 + b]++;
                        m &= m - 1;
                    }
                }
            }
            for (uint64_t k = 0; k < r.invCnt; k++) {
                const uint64_t s = invalidPos[k];
                if (s >= d.g0len)
                    continue;
                const uint64_t* mp = mask.data() + (size_t) s * maskWords;
                for (uint32_t w = 0; w < maskWords; w++) {
                    uint64_t m = mp[w];
                    while (m) {
                        const int b = __builtin_ctzll(m);
                        iRow[(size_t) w * 64 + b]++;
                        m &= m - 1;
                    }
                }
            }
        }
        for (size_t i = 0; i < N; i++) {
            string& ob = outBufs[i];
            for (size_t ri = 0; ri < recs.size(); ri++) {
                const QueryRec& r = recs[ri];
                ob.append(d.buf.data() + r.nameOff, r.nameLen);
                ob.push_back('\t');
                ob += to_string(diffMat[ri * N + i]);
                ob.push_back('\t');
                ob += to_string(inRangeSites[i] - invMat[ri * N + i]);
                ob.push_back('\n');
            }
            if (ob.size() >= (8u << 20)) {
                outs[i]->write(ob.data(), ob.size());
                ob.clear();
            }
        }
        recs.clear();
        posBuf.clear();
    };

    const uint8_t* p = (const uint8_t*) d.buf.data() + d.dataOff;
    const uint8_t* end = (const uint8_t*) d.buf.data() + d.buf.size();
    const uint8_t* bufStart = (const uint8_t*) d.buf.data();
    auto decodePosList = [&](uint64_t& off, uint64_t& cnt) {
        uint32_t n;
        if (!readLE<uint32_t>(p, end, n)) {
            fprintf(stderr, "ERROR: %s is corrupted (truncated position list).\n", diffBinFileName.c_str());
            exit(EXIT_FAILURE);
        }
        off = posBuf.size();
        cnt = n;
        posBuf.resize(posBuf.size() + n);
        uint64_t prev = 0, delta;
        for (uint32_t k = 0; k < n; k++) {
            if (!readVarint(p, end, delta)) {
                fprintf(stderr, "ERROR: %s is corrupted (truncated position list).\n", diffBinFileName.c_str());
                exit(EXIT_FAILURE);
            }
            prev += delta;
            posBuf[off + k] = prev;
        }
    };

    while (p < end) {
        uint16_t nameLen;
        if (!readLE<uint16_t>(p, end, nameLen) || (size_t) (end - p) < nameLen) {
            fprintf(stderr, "ERROR: %s is corrupted (truncated record).\n", diffBinFileName.c_str());
            exit(EXIT_FAILURE);
        }
        QueryRec r;
        r.nameOff = (uint64_t) (p - bufStart);
        r.nameLen = nameLen;
        p += nameLen;
        decodePosList(r.diffOff, r.diffCnt);
        decodePosList(r.invOff, r.invCnt);
        recs.push_back(r);
        if (posBuf.size() >= BLOCK_POS_LIMIT)
            flushBlock();
    }
    flushBlock();
    for (size_t i = 0; i < N; i++) {
        outs[i]->write(outBufs[i].data(), outBufs[i].size());
        if (outFiles[i].is_open())
            outFiles[i].close();
    }
}

} // namespace

void MBGC_Encoder::queryDiff(const string& diffBinFileName, const string& sitesFileName, const string& outFileName,
                             bool forceOverwrite) {
    vector<pair<string, string>> jobs;
    jobs.emplace_back(sitesFileName, outFileName);
    queryDiffBatch(diffBinFileName, jobs, forceOverwrite);
}

void MBGC_Encoder::queryDiffBatch(const string& diffBinFileName, const vector<pair<string, string>>& jobs,
                                  bool forceOverwrite) {
    LoadedDiff d;
    if (!loadDiffFile(diffBinFileName, d))
        exit(EXIT_FAILURE);

    int totalThreads = params->noOfThreadsLimited ? params->coderThreads : omp_get_num_procs();
    if (totalThreads < 1)
        totalThreads = 1;

    const size_t jobCount = jobs.size();
    if (jobCount == 1) {
        runDiffQueryJob(d, diffBinFileName, jobs[0], forceOverwrite, totalThreads);
        return;
    }

    // N >= 2: single diff load + single varint parse; one pass updates all
    // samples that own each position (union mask) - scales flat with N.
    runDiffQueryJobsUnion(d, diffBinFileName, jobs, forceOverwrite, totalThreads);
}

// ---------------------------------------------------------------------------
// preStrain representative-selection pipeline (single command, in-memory):
//  1) count sample diffs vs G0 at the given coords sites       -> ref_num
//  2) count diff_count of every strain in the diff DB at the same coords
//  3) cluster-aware selection: closest strain per cluster, top X globally
// ---------------------------------------------------------------------------

static string readFastaSeq(const string& path) {
    ifstream in(path, ios::binary);
    if (!in) {
        fprintf(stderr, "ERROR: cannot open FASTA file %s\n", path.c_str());
        exit(EXIT_FAILURE);
    }
    string seq, line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '>')
            continue;
        for (char& c : line)
            c = (char) toupper((unsigned char) c);
        line.erase(remove_if(line.begin(), line.end(),
                             [](char c) { return isspace((unsigned char) c) != 0; }), line.end());
        seq += line;
    }
    return seq;
}

static uint64_t computeRefNum(const string& g0Seq, const string& refFile, const string& coordsFile) {
    const string ref = readFastaSeq(refFile);
    if (g0Seq.size() != ref.size()) {
        fprintf(stderr, "ERROR: length mismatch - G0 %zu bp, really-ref.fa %zu bp (must be equal).\n",
                g0Seq.size(), ref.size());
        exit(EXIT_FAILURE);
    }
    ifstream cf(coordsFile);
    if (!cf) {
        fprintf(stderr, "ERROR: cannot open coords file %s\n", coordsFile.c_str());
        exit(EXIT_FAILURE);
    }
    uint64_t diff = 0;
    string line;
    while (getline(cf, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        size_t b = line.find_first_not_of(" \t");
        if (b == string::npos)
            continue;
        size_t e = line.find_last_not_of(" \t");
        if (line[b] == '#')
            continue;
        uint64_t pos = strtoull(line.substr(b, e - b + 1).c_str(), nullptr, 10); // 1-based
        if (pos == 0 || pos - 1 >= g0Seq.size())
            continue;
        const char a = g0Seq[pos - 1], bch = ref[pos - 1];
        const bool aN = (a == 'A' || a == 'C' || a == 'G' || a == 'T');
        const bool bN = (bch == 'A' || bch == 'C' || bch == 'G' || bch == 'T');
        if (aN && bN && a != bch)
            diff++;
    }
    return diff;
}

// Count, for every strain record in the diff stream, how many of its diff
// positions fall into the given sites (1-based, in-range). Runs in memory,
// one parse of the raw diff, parallel per-block counting.
static void queryDiffCounts(const LoadedDiff& d, const vector<uint64_t>& sites,
                            vector<string>& names, vector<uint64_t>& diffs) {
    uint64_t outOfRange = 0;
    for (uint64_t s : sites)
        if (s >= d.g0len)
            outOfRange++;
    if (outOfRange > 0)
        fprintf(stderr, "WARNING: %llu site(s) exceed the G0 sequence length (%llu bp) "
                        "and are excluded from the count.\n",
                (unsigned long long) outOfRange, (unsigned long long) d.g0len);

    vector<uint64_t> siteBits((size_t) ((d.g0len + 63) / 64), 0);
    for (uint64_t s : sites)
        if (s < d.g0len)
            siteBits[(size_t) (s >> 6)] |= 1ULL << (s & 63);

    struct QRec {
        uint64_t nameOff;
        uint32_t nameLen;
        uint64_t diffOff, diffCnt, invOff, invCnt;
    };
    const size_t BLOCK_POS_LIMIT = 1u << 24; // ~16M positions (~128 MB) per block
    vector<QRec> recs;
    vector<uint64_t> posBuf;

    auto flushBlock = [&](int threads) {
        if (recs.empty())
            return;
        vector<uint64_t> blockDiffs(recs.size(), 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 32)
        for (int64_t i = 0; i < (int64_t) recs.size(); i++) {
            const QRec& r = recs[i];
            const uint64_t* dp = posBuf.data() + r.diffOff;
            uint64_t c = 0;
            for (uint64_t k = 0; k < r.diffCnt; k++) {
                const uint64_t s = dp[k];
                if (s < d.g0len)
                    c += (siteBits[(size_t) (s >> 6)] >> (s & 63)) & 1ULL;
            }
            blockDiffs[i] = c;
        }
        for (size_t i = 0; i < recs.size(); i++) {
            names.emplace_back(d.buf.data() + recs[i].nameOff, recs[i].nameLen);
            diffs.push_back(blockDiffs[i]);
        }
        recs.clear();
        posBuf.clear();
    };

    const uint8_t* p = (const uint8_t*) d.buf.data() + d.dataOff;
    const uint8_t* end = (const uint8_t*) d.buf.data() + d.buf.size();
    const uint8_t* bufStart = (const uint8_t*) d.buf.data();

    auto decodePosList = [&](uint64_t& off, uint64_t& cnt) {
        uint32_t n;
        if (!readLE<uint32_t>(p, end, n)) {
            fprintf(stderr, "ERROR: diff is corrupted (truncated position list).\n");
            exit(EXIT_FAILURE);
        }
        off = posBuf.size();
        cnt = n;
        posBuf.resize(posBuf.size() + n);
        uint64_t prev = 0, delta;
        for (uint32_t k = 0; k < n; k++) {
            if (!readVarint(p, end, delta)) {
                fprintf(stderr, "ERROR: diff is corrupted (truncated position list).\n");
                exit(EXIT_FAILURE);
            }
            prev += delta;
            posBuf[off + k] = prev;
        }
    };

    const int threads = omp_get_num_procs();
    while (p < end) {
        uint16_t nameLen;
        if (!readLE<uint16_t>(p, end, nameLen) || (size_t) (end - p) < nameLen) {
            fprintf(stderr, "ERROR: diff is corrupted (truncated record).\n");
            exit(EXIT_FAILURE);
        }
        QRec r;
        r.nameOff = (uint64_t) (p - bufStart);
        r.nameLen = nameLen;
        p += nameLen;
        decodePosList(r.diffOff, r.diffCnt);
        decodePosList(r.invOff, r.invCnt);
        recs.push_back(r);
        if (posBuf.size() >= BLOCK_POS_LIMIT)
            flushBlock(threads);
    }
    flushBlock(threads);
}

static vector<string> selectRepresentatives(const vector<string>& names, const vector<uint64_t>& diffs,
                                            const string& clusterFile, uint64_t refNum, uint64_t pickX,
                                            string& repText) {
    unordered_map<string, uint32_t> clusterMap;
    ifstream cf(clusterFile);
    if (!cf) {
        fprintf(stderr, "ERROR: cannot open cluster file %s\n", clusterFile.c_str());
        exit(EXIT_FAILURE);
    }
    string line;
    while (getline(cf, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        const size_t comma = line.find(',');
        if (comma == string::npos)
            continue;
        clusterMap[line.substr(0, comma)] = (uint32_t) strtoull(line.c_str() + comma + 1, nullptr, 10);
    }

    struct Cand {
        string name;
        uint64_t absVal;
        uint32_t cid;
    };
    vector<Cand> cands;
    cands.reserve(names.size());
    for (size_t i = 0; i < names.size(); i++) {
        auto it = clusterMap.find(names[i]);
        if (it == clusterMap.end())
            continue;
        const uint64_t av = diffs[i] > refNum ? diffs[i] - refNum : refNum - diffs[i];
        cands.push_back({names[i], av, it->second});
    }
    sort(cands.begin(), cands.end(),
         [](const Cand& a, const Cand& b) { return a.absVal < b.absVal; });

    // per cluster: keep the closest strain (first seen after sorting), then top X globally
    unordered_set<uint32_t> used;
    vector<string> picked;
    picked.reserve(pickX);
    for (const Cand& c : cands) {
        if (used.count(c.cid))
            continue;
        used.insert(c.cid);
        repText += c.name + " " + to_string(c.absVal) + " " + to_string(c.cid) + "\n";
        picked.push_back(c.name);
        if (picked.size() >= pickX)
            break;
    }
    return picked;
}

vector<string> MBGC_Encoder::preStrainSelect(const string& g0File, const string& refFile,
                                             const string& coordsFile, const string& diffFile,
                                             const string& clusterFile, uint64_t pickX,
                                             uint64_t refNumArg, bool refNumProvided,
                                             string& repText) {
    const uint64_t refNum = refNumProvided ? refNumArg : computeRefNum(readFastaSeq(g0File), refFile, coordsFile);
    if (refNumProvided)
        fprintf(stderr, "sample diffs vs G0 at coords sites: %llu (provided via -n)\n", (unsigned long long) refNum);
    else
        fprintf(stderr, "sample diffs vs G0 at coords sites: %llu\n", (unsigned long long) refNum);

    vector<uint64_t> sites = loadSitesFile(coordsFile);
    LoadedDiff d;
    if (!loadDiffFile(diffFile, d))
        exit(EXIT_FAILURE);

    vector<string> names;
    vector<uint64_t> diffs;
    queryDiffCounts(d, sites, names, diffs);
    fprintf(stderr, "DB strains scored: %zu\n", names.size());

    return selectRepresentatives(names, diffs, clusterFile, refNum, pickX, repText);
}

string MBGC_Encoder::getDiffG0Name(const string& diffBinFileName) {
    LoadedDiff d;
    if (!loadDiffFile(diffBinFileName, d))
        exit(EXIT_FAILURE);
    return d.g0name;
}

/* batch preStrain: shared cluster map version of selectRepresentatives */
static vector<string> selectRepresentativesWithMap(const vector<string>& names, const vector<uint64_t>& diffs,
                                                   const unordered_map<string, uint32_t>& clusterMap,
                                                   uint64_t refNum, uint64_t pickX, string& repText) {
    struct Cand {
        string name;
        uint64_t absVal;
        uint32_t cid;
    };
    vector<Cand> cands;
    cands.reserve(names.size());
    for (size_t i = 0; i < names.size(); i++) {
        auto it = clusterMap.find(names[i]);
        if (it == clusterMap.end())
            continue;
        const uint64_t av = diffs[i] > refNum ? diffs[i] - refNum : refNum - diffs[i];
        cands.push_back({names[i], av, it->second});
    }
    sort(cands.begin(), cands.end(),
         [](const Cand& a, const Cand& b) { return a.absVal < b.absVal; });
    unordered_set<uint32_t> used;
    vector<string> picked;
    picked.reserve(pickX);
    for (const Cand& c : cands) {
        if (used.count(c.cid))
            continue;
        used.insert(c.cid);
        repText += c.name + " " + to_string(c.absVal) + " " + to_string(c.cid) + "\n";
        picked.push_back(c.name);
        if (picked.size() >= pickX)
            break;
    }
    return picked;
}

vector<MBGC_Encoder::PreSampleResult>
MBGC_Encoder::preStrainBatch(const string& g0File, const vector<PreSample>& samples,
                            const string& diffFile, const string& clusterFile,
                            uint64_t pickX) {
    const size_t N = samples.size();
    vector<PreSampleResult> results(N);

    // 1) load the diff DB exactly once (also gives g0len)
    auto t0 = std::chrono::steady_clock::now();
    LoadedDiff d;
    if (!loadDiffFile(diffFile, d))
        exit(EXIT_FAILURE);
    auto t1 = std::chrono::steady_clock::now();
    const size_t nwords = (size_t)((d.g0len + 63) / 64);

    // 2) per-sample ref_num + coords bitset (parallel across samples)
    // G0.fa read once, shared across all samples
    const string g0Seq = readFastaSeq(g0File);
    vector<vector<uint64_t>> sampleBits(N, vector<uint64_t>(nwords, 0));
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < N; i++) {
        results[i].name = samples[i].name;
        results[i].refNum = computeRefNum(g0Seq, samples[i].refFile, samples[i].coordsFile);
        vector<uint64_t> sites = loadSitesFile(samples[i].coordsFile);
        for (uint64_t s : sites)
            if (s < d.g0len)
                sampleBits[i][(size_t)(s >> 6)] |= 1ULL << (s & 63);
        fprintf(stderr, "sample %s: ref_num=%llu, sites=%zu\n",
                samples[i].name.c_str(), (unsigned long long)results[i].refNum, sites.size());
    }
    auto t2 = std::chrono::steady_clock::now();

    // 3) ONE streaming pass over the diff records: score every strain vs ALL samples
    vector<string> names;
    vector<vector<uint64_t>> counts(N);

    struct QRec {
        uint64_t nameOff;
        uint32_t nameLen;
        uint64_t diffOff, diffCnt, invOff, invCnt;
    };
    const size_t BLOCK_POS_LIMIT = 1u << 24;
    vector<QRec> recs;
    vector<uint64_t> posBuf;
    vector<vector<uint64_t>> blockDiffs(N);

    auto flushBlock = [&](int threads) {
        if (recs.empty())
            return;
        for (size_t s = 0; s < N; s++)
            blockDiffs[s].assign(recs.size(), 0);
#pragma omp parallel for num_threads(threads) schedule(dynamic, 32)
        for (int64_t i = 0; i < (int64_t) recs.size(); i++) {
            const QRec& r = recs[i];
            const uint64_t* dp = posBuf.data() + r.diffOff;
            for (uint64_t k = 0; k < r.diffCnt; k++) {
                const uint64_t s = dp[k];
                if (s >= d.g0len)
                    continue;
                const size_t w = (size_t)(s >> 6);
                const uint64_t bit = 1ULL << (s & 63);
                for (size_t si = 0; si < N; si++)
                    blockDiffs[si][i] += (sampleBits[si][w] & bit) ? 1ULL : 0ULL;
            }
        }
        for (size_t i = 0; i < recs.size(); i++)
            names.emplace_back(d.buf.data() + recs[i].nameOff, recs[i].nameLen);
        for (size_t s = 0; s < N; s++)
            counts[s].insert(counts[s].end(), blockDiffs[s].begin(), blockDiffs[s].end());
        recs.clear();
        posBuf.clear();
    };

    const uint8_t* p = (const uint8_t*) d.buf.data() + d.dataOff;
    const uint8_t* end = (const uint8_t*) d.buf.data() + d.buf.size();
    const uint8_t* bufStart = (const uint8_t*) d.buf.data();
    auto decodePosList = [&](uint64_t& off, uint64_t& cnt) {
        uint32_t n;
        if (!readLE<uint32_t>(p, end, n)) {
            fprintf(stderr, "ERROR: diff is corrupted (truncated position list).\n");
            exit(EXIT_FAILURE);
        }
        off = posBuf.size();
        cnt = n;
        posBuf.resize(posBuf.size() + n);
        uint64_t prev = 0, delta;
        for (uint32_t k = 0; k < n; k++) {
            if (!readVarint(p, end, delta)) {
                fprintf(stderr, "ERROR: diff is corrupted (truncated position list).\n");
                exit(EXIT_FAILURE);
            }
            prev += delta;
            posBuf[off + k] = prev;
        }
    };

    const int threads = omp_get_num_procs();
    while (p < end) {
        uint16_t nameLen;
        if (!readLE<uint16_t>(p, end, nameLen) || (size_t)(end - p) < nameLen) {
            fprintf(stderr, "ERROR: diff is corrupted (truncated record).\n");
            exit(EXIT_FAILURE);
        }
        QRec r;
        r.nameOff = (uint64_t)(p - bufStart);
        r.nameLen = nameLen;
        p += nameLen;
        decodePosList(r.diffOff, r.diffCnt);
        decodePosList(r.invOff, r.invCnt);
        recs.push_back(r);
        if (posBuf.size() >= BLOCK_POS_LIMIT)
            flushBlock(threads);
    }
    flushBlock(threads);
    auto t3 = std::chrono::steady_clock::now();
    fprintf(stderr, "batch scored %zu strains once for %zu samples\n", names.size(), N);
    fprintf(stderr, "TIMING: loadDiff=%.2fs refNum+bitset=%.2fs streamingScore=%.2fs\n",
            std::chrono::duration<double>(t1-t0).count(),
            std::chrono::duration<double>(t2-t1).count(),
            std::chrono::duration<double>(t3-t2).count());

    // 4) shared cluster map, per-sample selection
    unordered_map<string, uint32_t> clusterMap;
    {
        ifstream cf(clusterFile);
        if (!cf) {
            fprintf(stderr, "ERROR: cannot open cluster file %s\n", clusterFile.c_str());
            exit(EXIT_FAILURE);
        }
        string line;
        while (getline(cf, line)) {
            if (!line.empty() && line.back() == '\r')
                line.pop_back();
            const size_t comma = line.find(',');
            if (comma == string::npos)
                continue;
            clusterMap[line.substr(0, comma)] = (uint32_t) strtoull(line.c_str() + comma + 1, nullptr, 10);
        }
    }
    for (size_t s = 0; s < N; s++)
        results[s].reps = selectRepresentativesWithMap(names, counts[s], clusterMap,
                                                       results[s].refNum, pickX, results[s].repText);
    return results;
}
