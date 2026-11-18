#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <chrono>
#include <set>
#include <sys/stat.h>
#include <dirent.h>

#include <zstd.h>

#include "mbgccoder/MBGC_Encoder.h"
#include "mbgccoder/MBGC_Decoder.h"
#include "mbgccoder/MBGC_Params.h"

#include <omp.h>

#ifdef __MINGW32__
#include <fcntl.h>
#endif

#define RELEASE_DATE "2026-08-29"

using namespace std;


void printVersion(bool details) {
    string date = RELEASE_DATE;
    if (!details)
        date.resize(4);
    fprintf(stderr, "Multiple Bacteria Genome Compressor (MBGC) v%d.%d.%d (c) Tomasz Kowalski, Szymon Grabowski, %s\n",
            (int) MBGC_Params::MBGC_VERSION_MAJOR, (int) MBGC_Params::MBGC_VERSION_MINOR,
            (int) MBGC_Params::MBGC_VERSION_REVISION, date.c_str());
}

void printCommonOpts() {
    fprintf(stderr,
            "\t[-%c <noOfThreads>] set limit of used threads\n"
            "\t[-%c] ignore FASTA file paths (use only filenames)\n",
            MBGC_Params::NO_OF_THREADS_OPT, MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT);
    fprintf(stderr,
            "\t[-%c] redirect app output to stderr\n",
            MBGC_Params::APP_OUT_TO_CERR_OPT);
    fprintf(stderr,
            "\t[-%c] print full command help and exit\n"
            "\t[-%c] print version number and exit\n",
            MBGC_Params::HELP_OPT, MBGC_Params::VERSION_OPT);
}

#ifdef DEVELOPER_BUILD
void printCommonDevOpts() {
    fprintf(stderr, "[-%c] printout params\n", MBGC_Params::PRINTOUT_PARAMS_OPT);
    fprintf(stderr, "[-%c <logfilename>] printout params to tab-separated filelog\n", MBGC_Params::PARAMS_TO_LOGFILE_OPT);
}

void printCompressDeveloperOptions(bool repackCommand) {
    fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
    fprintf(stderr, "[-%c] brute parallel encoding mode (without producer-consumer approach)\n",
            MBGC_Params::BRUTE_PARALLEL_MODE_OPT);
    fprintf(stderr, "[-%c] stronger backend compression (slow)\n", MBGC_Params::ULTRA_STREAMS_COMPRESSION_OPT);
    fprintf(stderr, "[-%c] protein sequences compression profile\n", MBGC_Params::PROTEINS_PROFILE_OPT);
    fprintf(stderr, "[-%c <skipMargin>] (0 <= %c <= 255, default: %d)\n", MBGC_Params::SKIP_MARGIN_OPT, MBGC_Params::SKIP_MARGIN_OPT, MBGC_Params::DEFAULT_SKIP_MARGIN);
    fprintf(stderr, "[-%c] limit reference buffer size up to 2^32 bytes\n", MBGC_Params::LIMIT_32_BIT_REF_OPT);
    fprintf(stderr, "[-%c <referenceSlidingWindowFactor>] (2 <= %c <= 255, default: %d)\n",
            MBGC_Params::REF_SLIDING_WINDOW_FACTOR_OPT, MBGC_Params::REF_SLIDING_WINDOW_FACTOR_OPT, MBGC_Params::DEFAULT_REFERENCE_SLIDING_WINDOW_FACTOR);
    fprintf(stderr, "--- MISMATCHES ENCODING AND MATCHING SCORE OPTIONS ---\n");
    fprintf(stderr, "[-%c <gapDepthMismatchesEncoding>] (0 <= %c <= %d, default: %d, 0 - disable)\n"
                    "[-%c <mismatchesLimitInPercent>] (0 < %c < 100, default: %d)\n"
                    "[-%c <initialScoreInPercent>] (0 <= %c <= 100, default: %d)\n",
            MBGC_Params::GAP_DEPTH_MISMATCHES_ENCODING_OPT, MBGC_Params::GAP_DEPTH_MISMATCHES_ENCODING_OPT, MBGC_Params::MAX_GAP_DEPTH, MBGC_Params::DEFAULT_GAP_DEPTH_MISMATCHES_ENCODING,
            MBGC_Params::MISMATCHES_LIMIT_IN_PERCENT_OPT, MBGC_Params::MISMATCHES_LIMIT_IN_PERCENT_OPT, MBGC_Params::DEFAULT_MISMATCHES_LIMIT_IN_PERCENT,
            MBGC_Params::MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT, MBGC_Params::MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT, MBGC_Params::DEFAULT_MISMATCHES_INITIAL_SCORE_IN_PERCENT);
    fprintf(stderr, "--- EXPERIMENTAL REF EXTENSION STRATEGIES ---\n");
    fprintf(stderr, "[-%c] combined ref extension strategy\n", MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_OPT);
    fprintf(stderr, "[-%c] literals ref extension strategy (experimental - only for encoding)\n"
                    "[-%c <literalContextLength>] (0 <= %c, default: %ld)\n"
                    "[-%c <minimalLiteralLength>] (1 <= %c, default: DISABLED)\n",
            MBGC_Params::ONLY_LITERAL_REF_EXTENSION_STRATEGY_OPT,
            MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT_OPT, MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT_OPT, MBGC_Params::DEFAULT_REF_LITERAL_BEFORE_AFTER_EXT,
            MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_OPT, MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_OPT);
    fprintf(stderr, "[-%c] split contigs into blocks (experimental - only for encoding)\n", MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_OPT);
    fprintf(stderr, "[-%c] dynamic unmatched fraction ref extension strategy\n", MBGC_Params::DYNAMIC_REF_EXTENSION_FACTOR_OPT);
    fprintf(stderr, "--- DISSIMILAR CONTIGS TWEAK OPTIONS ---\n");
    fprintf(stderr, "[-%c <allowedTargetsOutrun>] (default: %d)\n",
            MBGC_Params::ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT, MBGC_Params::DEFAULT_ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS);
    fprintf(stderr, "[-%c <minimalDissimilarContigLength>] (default: %ld)\n",
            MBGC_Params::MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT, MBGC_Params::DEFAULT_MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS);
    fprintf(stderr, "[-%c <dissimilarContigUnmatchedFractionTweakFactor>] (default: %d)\n",
            MBGC_Params::UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT, MBGC_Params::DEFAULT_UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS);
    fprintf(stderr, "--- INPUT/OUTPUT OPTIONS ---\n");
    if (!repackCommand)
        fprintf(stderr, "[-%c] interleave input files order\n", MBGC_Params::INTERLEAVE_FILE_ORDER_OPT);
    printCommonDevOpts();
}
#endif

void printEncodeUsage(string base_cmd_name, bool details, bool repackMode) {
    printVersion(details);
    if (repackMode) {
        fprintf(stderr, "\nUsage for repacking part of archive:\n\t%s [-%c <compressionMode>]"
                        " [-%c <pattern>] [-%c <patternsFile>]\n\t\t<archiveFile> <outputArchiveFile>\n", base_cmd_name.c_str(),
                        MBGC_Params::COMPRESSION_MODE_OPT,
                        MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
        if (details) {
            fprintf(stderr,
                    "\n<archiveFile> mbgc archive filename for repacking\n"
                    "\tfor standard input, set <archiveFile> to %s\n"
                    "<patternsFile> name of text file with list of patterns (in separate lines)\n"
                    "\texcludes files not matching any pattern (does not invalidate -%c option)\n"
                    "<outputArchiveFile> name of repacked mbgc archive\n"
                    "\tfor standard output, set <outputArchiveFile> to %s\n",
                    MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::STANDARD_IO_POSIX_ALIAS);
        }
    } else {
        fprintf(stderr, "\nUsage for multiple file compression (list of files given as input):\n\t%s [-%c <compressionMode>]"
                        " [-%c] <fastaFileList> <archiveFile>\n", base_cmd_name.c_str(),
                MBGC_Params::COMPRESSION_MODE_OPT, MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT);
        fprintf(stderr, "Usage for single file compression:\n\t%s [-%c <compressionMode>]"
                        " [-%c] -%c <inputFastaFile> <archiveFile>\n", base_cmd_name.c_str(),
                MBGC_Params::COMPRESSION_MODE_OPT, MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT, MBGC_Params::INPUT_FILE_OPT);

        if (details) {
            fprintf(stderr,
                    "\n<fastaFileList> name of text file with a list of FASTA files (raw or gz)\n"
                    "\t(given in separate lines) for compression\n"
                    "<inputFastaFile> name of a FASTA file (raw or gz) for compression\n"
                    "\tfor standard input, set <inputFastaFile> to %s\n"
                    "<archiveFile> mbgc archive filename\n"
                    "\tfor standard output in compression, set <archiveFile> to %s\n",
                    MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS);
        }
    }
    fprintf(stderr, "\nBasic options:\n\t[-%c <compressionMode>] "
                    "(speed: 0; default: 1; repo: 2; max: 3)\n"
                    "\t[-%c] allow lossy compression\n"
                    "\t[-%c] overwrite an existing output file\n",
                    MBGC_Params::COMPRESSION_MODE_OPT, MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT,
                    MBGC_Params::FORCE_OVERWRITE_OPT);

    if (repackMode) {
        fprintf(stderr,
                "\t[-%c <basesPerRow>] custom format of repacked DNA (0 - unlimited)\n"
                "\t[-%c <pattern>] exclude files with names not containing pattern\n"
                "\t[-%c <patternsFile>] exclude files not matching any pattern\n",
                MBGC_Params::DNA_LINE_LENGTH_OPT, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
        fprintf(stderr,
            "\t[-%c] disable parallel coding (does not apply to I/O and backend)\n",
            MBGC_Params::SEQ_CODING_MODE_OPT);
    } else {
        fprintf(stderr,
                "\t[-%c] disable parallel matching (does not apply to I/O and backend)\n",
                MBGC_Params::SEQ_CODING_MODE_OPT);
    }
    printCommonOpts();
    if (details) {
        fprintf(stderr, "\nCompression modes description:\n"
                        "\t(0) speed - for speed (fastest compression and decompression)\n"
                        "\t(1) default - regular mode (good ratio, fast)\n"
                        "\t(2) repo - for public repositories (better ratio, good speed)\n"
                        "\t(3) max - for long-term storage (best ratio, memory-frugal)\n");
    }
}

void printEncodeAdvancedDetails() {
    fprintf(stderr, "\n------------------ ADVANCED OPTIONS ----------------\n");
    fprintf(stderr, "[-%c <noOfWorkers>] worker threads used by matcher (default: %d)\n",
                    MBGC_Params::MATCHER_NO_OF_THREADS_OPT, MBGC_Params::DEFAULT_MATCHER_NO_OF_THREADS);
    fprintf(stderr, "[-%c] converts bases to uppercase\n",
            MBGC_Params::UPPERCASE_DNA_OPT);
    fprintf(stderr, "[-%c <matchingKmerLength>] (%d <= %c <= 40, default: %d)\n"
                    "[-%c <unmatchedFractionFactor>] (1 <= %c <= 255, default: %d)\n",
            MBGC_Params::KMER_LENGTH_OPT,
            MBGC_Params::MIN_MATCH_LENGTH, MBGC_Params::KMER_LENGTH_OPT, MBGC_Params::DEFAULT_KMER_LENGTH,
            MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT, MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT, MBGC_Params::DEFAULT_UNMATCHED_LENGTH_FACTOR);
    fprintf(stderr, "[-%c <unmatchedFractionRCFactor>] (0 <= %c <= 255, default: %d,\n\t\t0 disables rc-matches in reference)\n",
            MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT, MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT, MBGC_Params::DEFAULT_UNMATCHED_LENGTH_RC_FACTOR);
    fprintf(stderr, "[-%c <referenceSamplingStep>] (%c > 0, default: %d)\n"
                    "[-%c <referenceFactorBinaryOrder>] (0 <= %c <= 12, auto-adjusted by default)\n"
                    "[-%c] disable circular reference buffer\n",
            MBGC_Params::REF_SAMPLING_STEP_OPT, MBGC_Params::REF_SAMPLING_STEP_OPT, MBGC_Params::DEFAULT_MODE_REF_SAMPLING,
            MBGC_Params::REF_FACTOR_BIN_ORDER_OPT, MBGC_Params::REF_FACTOR_BIN_ORDER_OPT, MBGC_Params::CIRCULAR_REFERENCE_BUFFER_OPT);
    fprintf(stderr,
            "[-%c <gapDepthOffsetEncoding>] (0 <= %c <= %d, default: %d, 0 - disable)\n"
            "[-%c <gapBreakingMatchMinimalLength>] (%c > 0, default: %d, 0 - disable)\n"
            "[-%c <maxConsecutiveMismatches>] (0 <= %c <= 255, default: %d; 0 - disable;\n\t\t1 - single mismatch mode)\n"
            "[-%c] disable encoding mismatches with exclusion\n"
            "[-%c <rcMatchMinimalLength>] (%c >= 24, default: %d, 0 - disable)\n",
            MBGC_Params::GAP_DEPTH_OFFSET_ENCODING_OPT, MBGC_Params::GAP_DEPTH_OFFSET_ENCODING_OPT, MBGC_Params::MAX_GAP_DEPTH, MBGC_Params::DEFAULT_GAP_DEPTH_OFFSET_ENCODING,
            MBGC_Params::GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT, MBGC_Params::GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT, MBGC_Params::DEFAULT_GAP_BREAKING_MATCH_MIN_LENGTH,
            MBGC_Params::MAX_CONSECUTIVE_MISMATCHES_OPT, MBGC_Params::MAX_CONSECUTIVE_MISMATCHES_OPT, MBGC_Params::DEFAULT_MAX_CONSECUTIVE_MISMATCHES,
            MBGC_Params::DISABLE_MISMATCHES_WITH_EXCLUSION_OPT,
            MBGC_Params::RC_MATCH_MINIMAL_LENGTH_OPT, MBGC_Params::RC_MATCH_MINIMAL_LENGTH_OPT, 0);
}

int encodeCommand(std::string base_cmd_name, int argc, char *argv[]) {
    base_cmd_name += string(" ") + string(argv[0]);
    char cmdId = argv[0][0];
    argv[0] = (char *) base_cmd_name.c_str();
    MBGC_Params params;
    params.repackCommand = cmdId == MBGC_Params::REPACK_CMD;
    const string SHORT_OPTS = params.repackCommand ? MBGC_Params::REPACK_SHORT_OPTS : MBGC_Params::COMPRESS_SHORT_OPTS;
    if (argc == 1) {
        printEncodeUsage(base_cmd_name, false, params.repackCommand);
        fprintf(stderr, "\n");
        return EXIT_SUCCESS;
    }
    int opt; // current option
    bool cliError = false;
    bool printParamsFlag = false;
    params.coderMode = MBGC_Params::DEFAULT_CODER_MODE;
    MBGC_Params repackInputParams;

    while ((opt = getopt(argc, argv, SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                params.forceOverwrite = true;
                break;
            case MBGC_Params::ALLOW_LOSSY_COMPRESSION_OPT:
                params.allowLossyParsing = true;
                break;
            case MBGC_Params::NO_OF_THREADS_OPT:
                params.setThreadsLimit(atoi(optarg));
                if (params.repackCommand)
                    repackInputParams.setThreadsLimit(atoi(optarg));
                break;
            case MBGC_Params::MATCHER_NO_OF_THREADS_OPT:
                params.setMatcherThreads(atoi(optarg));
                break;
            case MBGC_Params::COMPRESSION_MODE_OPT:
                params.setCompressionMode(atoi(optarg));
                break;
            case MBGC_Params::INPUT_FILE_OPT:
                params.setInputFileName(string(optarg));
                break;
            case MBGC_Params::SEQ_CODING_MODE_OPT:
                params.setSequentialMatchingMode();
                if (params.repackCommand)
                    repackInputParams.disableLazyDecompression = true;
                break;
            case MBGC_Params::KMER_LENGTH_OPT:
                params.setKmerLength(atoi(optarg));
                break;
            case MBGC_Params::UNMATCHED_FRACTION_FACTOR_OPT:
                params.setUnmatchedFractionFactor(atoi(optarg));
                break;
            case MBGC_Params::REF_SAMPLING_STEP_OPT:
                params.setReferenceSamplingStep(atoi(optarg));
                break;
            case MBGC_Params::REF_FACTOR_BIN_ORDER_OPT:
                params.setReferenceFactorBinaryOrder(atoi(optarg));
                break;
            case MBGC_Params::CIRCULAR_REFERENCE_BUFFER_OPT:
                params.disableCircularReference();
                break;
            case MBGC_Params::UNMATCHED_FRACTION_RC_FACTOR_OPT:
                params.setUnmatchedFractionRCFactor(atoi(optarg));
                break;
            case MBGC_Params::GAP_DEPTH_OFFSET_ENCODING_OPT:
                params.setGapDepthOffsetEncoding(atoi(optarg));
                break;
            case MBGC_Params::GAP_BREAKING_MATCH_MINIMAL_LENGTH_OPT:
                params.setGapBreakingMatchMinimalLength(atoi(optarg));
                break;
            case MBGC_Params::MAX_CONSECUTIVE_MISMATCHES_OPT:
                params.setMaxConsecutiveMismatches(atoi(optarg));
                break;
            case MBGC_Params::DISABLE_MISMATCHES_WITH_EXCLUSION_OPT:
                params.disableMismatchesEncodingWithExclusion();
                break;
            case MBGC_Params::RC_MATCH_MINIMAL_LENGTH_OPT:
                params.setRCMatchMinimalLength(atoi(optarg));
                break;
            case MBGC_Params::FILTER_PATTERN_OPT:
                repackInputParams.setFilterPattern(string(optarg));
                break;
            case MBGC_Params::PATTERNS_FILE_OPT:
                repackInputParams.setInputFileName(string(optarg));
                break;
            case MBGC_Params::DNA_LINE_LENGTH_OPT:
                repackInputParams.setDNALineLength(atoll(optarg));
                break;
            case MBGC_Params::UPPERCASE_DNA_OPT:
                params.uppercaseDNA = true;
                break;
#ifdef DEVELOPER_BUILD
            case MBGC_Params::BRUTE_PARALLEL_MODE_OPT:
                params.setBruteParallelMode();
                break;
            case MBGC_Params::ULTRA_STREAMS_COMPRESSION_OPT:
                params.setUltraStreamsCompression();
                break;
            case MBGC_Params::SKIP_MARGIN_OPT:
                params.setSkipMargin(atoi(optarg));
                break;
            case MBGC_Params::LIMIT_32_BIT_REF_OPT:
                params.limit32bitReference();
                break;
            case MBGC_Params::REF_SLIDING_WINDOW_FACTOR_OPT:
                params.setReferenceSlidingWindowFactor(atoi(optarg));
                break;
            case MBGC_Params::GAP_DEPTH_MISMATCHES_ENCODING_OPT:
                params.setGapDepthMismatchesEncoding(atoi(optarg));
                break;
            case MBGC_Params::MISMATCHES_LIMIT_IN_PERCENT_OPT:
                params.setMismatchesLimitInPercent(atoi(optarg));
                break;
            case MBGC_Params::MISMATCHES_INITIAL_SCORE_IN_PERCENT_OPT:
                params.setMismatchesInitialScoreInPercent(atoi(optarg));
                break;
            case MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_OPT:
                params.refExtensionStrategy |= MBGC_Params::COMBINED_REF_EXTENSION_STRATEGY_MASK;
                break;
            case MBGC_Params::ONLY_LITERAL_REF_EXTENSION_STRATEGY_OPT:
                params.refExtensionStrategy |= MBGC_Params::ONLY_LITERAL_REF_EXTENSION_STRATEGY_MASK;
                break;
            case MBGC_Params::REF_LITERAL_BEFORE_AFTER_EXT_OPT:
                params.refLiteralBeforeAfterExt = atoi(optarg);
                break;
            case MBGC_Params::REF_LITERAL_MINIMAL_LENGTH_OPT:
                params.refLiteralMinimalLengthExt = atoi(optarg);
                break;
            case MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_OPT:
                params.refExtensionStrategy |= MBGC_Params::BLOCK_REF_EXTENSION_STRATEGY_MASK;
                break;
            case MBGC_Params::DYNAMIC_REF_EXTENSION_FACTOR_OPT:
                params.refExtensionStrategy |= MBGC_Params::DYNAMIC_REF_EXT_FACTOR_MASK;
                break;
            case MBGC_Params::INTERLEAVE_FILE_ORDER_OPT:
                params.setInterleaveFileOrder();
                break;
            case MBGC_Params::ALLOWED_TARGETS_OUTRUN_FOR_DISSIMILAR_CONTIGS_OPT:
                params.allowedTargetsOutrunForDissimilarContigs = atoi(optarg);
                break;
            case MBGC_Params::MINIMAL_LENGTH_FOR_DISSIMILAR_CONTIGS_TWEAK_OPT:
                params.minimalLengthForDissimilarContigs = atoi(optarg);
                break;
            case MBGC_Params::UNMATCHED_FRACTION_TWEAK_FACTOR_FOR_DISSIMILAR_CONTIGS_OPT:
                params.unmatchedFractionFactorTweakForDissimilarContigs = atoi(optarg);
                break;
            case MBGC_Params::PROTEINS_PROFILE_OPT:
                params.setProteinsCompressionProfile();
                break;
            case MBGC_Params::CONCAT_HEADERS_AND_SEQS_OPT:
                params.concatHeadersAndSequencesMode = true;
                break;
            case MBGC_Params::PRINTOUT_PARAMS_OPT:
                printParamsFlag = true;
                break;
            case MBGC_Params::PARAMS_TO_LOGFILE_OPT:
                params.enableLogFile(optarg);
                break;
#endif
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::APP_OUT_TO_CERR_OPT:
                params.setAppOutToCerr();
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printEncodeUsage(base_cmd_name, true, params.repackCommand);
                printEncodeAdvancedDetails();
#ifdef DEVELOPER_BUILD
                printCompressDeveloperOptions(params.repackCommand);
#endif
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n\n");
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
    if (params.repackCommand) {
        if (optind > (argc - 2) || optind < (argc - 2)) {
            fprintf(stderr, "%s: For 'repack' command expected 2 arguments after options "
                            "(found %d)\n", base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
    } else {
        if (params.inputFileName.empty() && (optind > (argc - 2) || optind < (argc - 2))) {
            fprintf(stderr, "%s: For 'compress' command using files list as input expected 2 arguments after options "
                            "(found %d)\n", base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
        if (!params.inputFileName.empty() && (optind > (argc - 1) || optind < (argc - 1))) {
            fprintf(stderr, "%s: For 'compress' command given single file as input expected 1 argument after options "
                            "(found %d)\n", base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

    if (!params.repackCommand && params.inputFileName.empty())
        params.setSeqListFileName(argv[optind++]);
    if (params.repackCommand)
        repackInputParams.setInArchiveFileName(argv[optind++]);
    params.setOutArchiveFileName(argv[optind++]);
    if (params.appOutToCerr) {
        PgHelpers::appout = &std::cerr;
#ifdef DEVELOPER_BUILD
        PgHelpers::devout = &std::cerr;
#endif
    } else if (params.outArchiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
        PgHelpers::devout = &null_stream;
    }
    MBGC_Encoder encoder(&params);
    if (params.repackCommand)
        encoder.repack(repackInputParams);
    else
        encoder.encode();
    *PgHelpers::appout << "encoding time - " << PgHelpers::time_millis() << " [ms]" << endl;

    if (printParamsFlag)
        params.printout();
    return EXIT_SUCCESS;
}

#ifdef DEVELOPER_BUILD
void printValidateUsage(string base_cmd_name, bool details) {
    printVersion(details);
    fprintf(stderr, "\nUsage for validation:\n\t%s [-%c <noOfThreads>] [-%c <dnaLineLength>]"
                    " <archiveFile>\n", base_cmd_name.c_str(),
            MBGC_Params::NO_OF_THREADS_OPT, MBGC_Params::DNA_LINE_LENGTH_OPT);
    fprintf(stderr, "Usage for partial validation:\n\t%s [-%c <pattern>] [-%c <patternsFile>]"
                    " <archiveFile>\n", base_cmd_name.c_str(), MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::INPUT_FILE_OPT);

    if (details) {
        fprintf(stderr,
                "\n<archiveFile> mbgc archive filename\n"
                "\tfor standard input in validation, set <archiveFile> to %s\n"
                "<patternsFile> name of text file with a list of patterns (in separate lines)\n"
                "\tfiles matching one of the patterns are validated obligatory (does not invalidate -%c option)\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::FILTER_PATTERN_OPT);
    }
    fprintf(stderr, "\nBasic options:\n"
                    "\t[-%c <basesPerRow>] custom format of decoded sequences (0 - unlimited)\n"
                    "\t[-%c <pattern>] validate files with names containing pattern\n"
                    "\t[-%c <patternsFile>] validate files with names containing any pattern\n",
            MBGC_Params::DNA_LINE_LENGTH_OPT, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
    fprintf(stderr, "\t[-%c] converts bases to uppercase\n",
            MBGC_Params::UPPERCASE_DNA_OPT);
    printCommonOpts();
}

void printDecompressDeveloperOptions(bool isListCmd, bool isValidationCmd) {
    fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
    if (isValidationCmd) {
        fprintf(stderr, "[-%c] dump streams to files\n", MBGC_Params::DUMP_STREAMS_OPT);
        fprintf(stderr, "[-%c] skip actual validation (measure decoding only)\n", MBGC_Params::SKIP_ACTUAL_VALIDATION_OPT);
    }
    if (!isListCmd && !isValidationCmd)
        fprintf(stderr, "[-%c] concatenate sequences and headers in each file\n", MBGC_Params::CONCAT_HEADERS_AND_SEQS_OPT);
    printCommonDevOpts();
}
#endif

void printDecodingUsage(string base_cmd_name, bool details, bool isListCmd) {
#ifdef DEVELOPER_BUILD
    if (base_cmd_name.back() == MBGC_Params::VALIDATE_CMD) {
        printValidateUsage(base_cmd_name, false);
        return;
    }
#endif
    printVersion(details);
    if (isListCmd) {
        fprintf(stderr, "\nUsage for partial file listing:\n\t%s [-%c] [-%c <pattern>] [-%c <patternsFile>]"
                        " <archiveFile>\n", base_cmd_name.c_str(),
                        MBGC_Params::LIST_HEADERS_OPT, MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
    } else {
        fprintf(stderr, "\nUsage for decompression:\n\t%s [-%c <gzLevel>]"
                        " <archiveFile> [<outputPath>]\n", base_cmd_name.c_str(), MBGC_Params::GZ_DECOMPRESSION_OPT);
        fprintf(stderr, "Usage for partial decompression (list of patterns given as input):\n\t%s [-%c <patternsFile>]"
                        " <archiveFile> [<outputPath>]\n", base_cmd_name.c_str(), MBGC_Params::PATTERNS_FILE_OPT);
    }

    if (details) {
        fprintf(stderr,
                "\n<archiveFile> mbgc archive filename\n"
                "\tfor standard input in decompression, set <archiveFile> to %s\n"
                "<patternsFile> name of text file with list of patterns (in separate lines)\n"
                "\texcludes files not matching any pattern (does not invalidate -%c option)\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::FILTER_PATTERN_OPT);
        if (!isListCmd)
            fprintf(stderr,
                    "<outputPath> extraction target path root (current directory by default)\n"
                    "\tfor standard output, set <outputPath> to %s (all files are concatenated)\n",
                    MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }
    fprintf(stderr, "\nBasic options:\n");
    if (isListCmd) {
        fprintf(stderr,
                "\t[-%c <pattern>] exclude files with names not containing pattern\n"
                "\t[-%c <patternsFile>] exclude files not matching any pattern\n"
                "\t[-%c] list sequence headers (using convention: \">sequencename>filename\")\n",
                MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT, MBGC_Params::LIST_HEADERS_OPT);
    } else {
        fprintf(stderr,
                "\t[-%c] overwrite an existing output files\n"
                "\t[-%c <gzLevel>] extract FASTA files to gz archives\n\t\t(compression level: 1 <= %c <= 12, recommended: %d)\n"
                "\t[-%c <basesPerRow>] custom format of decoded sequences (0 - unlimited)\n"
                "\t[-%c <pattern>] exclude files with names not containing pattern\n"
                "\t[-%c <patternsFile>] exclude files not matching any pattern\n",
                MBGC_Params::FORCE_OVERWRITE_OPT, MBGC_Params::GZ_DECOMPRESSION_OPT, MBGC_Params::GZ_DECOMPRESSION_OPT,
                MBGC_Params::RECOMMENDED_GZ_COMPRESSION_LEVEL, MBGC_Params::DNA_LINE_LENGTH_OPT,
                MBGC_Params::FILTER_PATTERN_OPT, MBGC_Params::PATTERNS_FILE_OPT);
        fprintf(stderr,
            "\t[-%c] disable parallel decoding (does not apply to I/O and backend)\n",
            MBGC_Params::SEQ_CODING_MODE_OPT);
    }
    printCommonOpts();
}

int decodeCommands(std::string base_cmd_name, int argc, char *argv[]) {
    string defaultCommandWarning;
    if (strlen(argv[0]) > 1)
        defaultCommandWarning = string("\nNote: selected default command '") + MBGC_Params::INFO_CMD + "'. Please use '" + base_cmd_name + "' to list all commands.\n";
    char cmdId = strlen(argv[0]) > 1 ? MBGC_Params::INFO_CMD : argv[0][0];
    base_cmd_name += string(" ") + cmdId;
    argv[0] = (char *) base_cmd_name.c_str();
    MBGC_Params params;
    params.infoCommand = cmdId == MBGC_Params::INFO_CMD;
#ifndef DEVELOPER_BUILD
    const string SHORT_OPTS = params.infoCommand ? MBGC_Params::INFO_SHORT_OPTS : MBGC_Params::DECOMPRESS_SHORT_OPTS;
#else
    if (params.infoCommand)
        PgHelpers::devout = &std::cerr;
    bool isValidateCmd = cmdId == MBGC_Params::VALIDATE_CMD;
    params.validationMode = isValidateCmd;
    const string SHORT_OPTS = isValidateCmd ? MBGC_Params::VALIDATION_SHORT_OPTS :
                              params.infoCommand ? MBGC_Params::INFO_SHORT_OPTS : MBGC_Params::DECOMPRESS_SHORT_OPTS;
#endif
    if (argc == 1) {
        printDecodingUsage(base_cmd_name, false, params.infoCommand);
        fprintf(stderr, "%s\n", defaultCommandWarning.c_str());
        return EXIT_SUCCESS;
    }
    int opt; // current option
    bool cliError = false;
    bool printParamsFlag = false;

    while ((opt = getopt(argc, argv, SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                params.forceOverwrite = true;
                break;
            case MBGC_Params::NO_OF_THREADS_OPT:
                params.setThreadsLimit(atoi(optarg));
                break;
            case MBGC_Params::FILTER_PATTERN_OPT:
                params.setFilterPattern(string(optarg));
                break;
            case MBGC_Params::PATTERNS_FILE_OPT:
                params.setInputFileName(string(optarg));
                break;
            case MBGC_Params::GZ_DECOMPRESSION_OPT:
                params.setGzDecompressionMode(atoi(optarg));
                break;
            case MBGC_Params::DNA_LINE_LENGTH_OPT:
                params.setDNALineLength(atoll(optarg));
                break;
            case MBGC_Params::LIST_HEADERS_OPT:
                params.listHeadersMode = true;
                break;
            case MBGC_Params::SEQ_CODING_MODE_OPT:
                params.disableLazyDecompression = true;
                break;
#ifdef DEVELOPER_BUILD
            case MBGC_Params::UPPERCASE_DNA_OPT:
                params.uppercaseDNA = true;
                break;
            case MBGC_Params::CONCAT_HEADERS_AND_SEQS_OPT:
                params.concatHeadersAndSequencesMode = true;
                break;
            case MBGC_Params::DUMP_STREAMS_OPT:
                dump_after_decompression = true;
                break;
            case MBGC_Params::SKIP_ACTUAL_VALIDATION_OPT:
                params.skipActualValidation = true;
                break;
            case MBGC_Params::PRINTOUT_PARAMS_OPT:
                printParamsFlag = true;
                break;
            case MBGC_Params::PARAMS_TO_LOGFILE_OPT:
                params.enableLogFile(optarg);
                break;
#endif
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::APP_OUT_TO_CERR_OPT:
                params.setAppOutToCerr();
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printDecodingUsage(base_cmd_name, true, params.infoCommand);
#ifdef DEVELOPER_BUILD
                printDecompressDeveloperOptions(params.infoCommand, isValidateCmd);
#endif
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n");
                fprintf(stderr, "%s\n", defaultCommandWarning.c_str());
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
#ifdef DEVELOPER_BUILD
    if (isValidateCmd && (optind < (argc - 1) || optind > (argc - 1))) {
        fprintf(stderr, "%s: For 'validation' command expected 1 argument after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    } else
#endif
    if (params.infoCommand && (optind < (argc - 1) || optind > (argc - 1))) {
        fprintf(stderr, "%s: For 'info' command expected 1 argument after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    } else if (optind < (argc - 2) || optind > (argc - 1)) {
        fprintf(stderr, "%s: For 'decompress' command expected 1 or 2 arguments after options (found %d)\n", argv[0],
                argc - optind);
        cliError = true;
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information%s\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT,
                defaultCommandWarning.c_str());
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

#ifdef DEVELOPER_BUILD
    dump_after_decompression_prefix = string(argv[optind]) + "_dump_";
#endif
    params.setInArchiveFileName(argv[optind++]);
    params.setOutputPath(argc > optind ? argv[optind++] : "");
    if (params.appOutToCerr) {
        PgHelpers::appout = &std::cerr;
#ifdef DEVELOPER_BUILD
        PgHelpers::devout = &std::cerr;
#endif
    } else  if (params.outputPath == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
        PgHelpers::devout = &null_stream;
    }
    *PgHelpers::logout << "            \t";
    MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&params, defaultCommandWarning);
    decoder->decode();
    delete(decoder);
    if (params.infoCommand)
        *PgHelpers::devout << "listing" << " time - " << PgHelpers::time_millis() << " [ms]" << endl;
    else
        *PgHelpers::appout << "decoding" << " time - " << PgHelpers::time_millis() << " [ms]" << endl;

    if (printParamsFlag) {
        params.printout();
    }
    return EXIT_SUCCESS;
}

#ifdef DEVELOPER_BUILD
void printAppendDeveloperOptions() {
    fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
    fprintf(stderr, "[-%c <matcherNumberOfThreads>] (default: %d)\n",
            MBGC_Params::MATCHER_NO_OF_THREADS_OPT, MBGC_Params::DEFAULT_MATCHER_NO_OF_THREADS);
    printCommonDevOpts();
}
#endif

void printAppendUsage(string base_cmd_name, bool details) {
    printVersion(details);
    fprintf(stderr, "\nUsage for appending multiple FASTA archive (list of files given as input):\n\t%s"
                    " <fastaFileList> <archiveFile> [<outputArchiveFile>]\n", base_cmd_name.c_str());
    fprintf(stderr, "Usage for appending single FASTA file archive:\n\t%s"
                    " -%c <inputFastaFile> <archiveFile> [<outputArchiveFile>]\n", base_cmd_name.c_str(), MBGC_Params::INPUT_FILE_OPT);

    if (details) {
        fprintf(stderr,
                "\n<fastaFileList> name of text file with a list of FASTA files (raw or gz)\n"
                "\t(given in separate lines) for appending to the <archiveFile>\n"
                "<inputFastaFile> name of a FASTA file (raw or gz) appended to the <archiveFile>\n"
                "\tfor standard input, set <inputFastaFile> to %s\n"
                "<archiveFile> mbgc archive filename to be appended\n"
                "\tfor standard input (and output if <outputArchiveFile> is not defined)\n"
                "\tset <archiveFile> to %s\n"
                "<outputArchiveFile> if defined, <archiveFile> archive remains unchanged\n"
                "\tand new archive is created\n"
                "\tfor standard output, set <outputArchiveFile> to %s\n",
                MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS, MBGC_Params::STANDARD_IO_POSIX_ALIAS);
    }
    fprintf(stderr, "\nBasic options:\n");
    fprintf(stderr, "\t[-%c] overwrite an existing output file\n", MBGC_Params::FORCE_OVERWRITE_OPT);
    fprintf(stderr,
            "\t[-%c] disable parallel coding (does not apply to I/O and backend)\n",
            MBGC_Params::SEQ_CODING_MODE_OPT);
    printCommonOpts();
}

int appendCommand(std::string base_cmd_name, int argc, char *argv[]) {
    base_cmd_name += string(" ") + string(argv[0]);
    argv[0] = (char *) base_cmd_name.c_str();
    if (argc == 1) {
        printAppendUsage(base_cmd_name, false);
        fprintf(stderr, "\n");
        return EXIT_SUCCESS;
    }
    int opt; // current option
    bool cliError = false;
    bool printParamsFlag = false;
    MBGC_Params params, newParams;
    newParams.coderThreads = -1;
    newParams.matcherWorkingThreads = -1;
    params.coderMode = MBGC_Params::DEFAULT_CODER_MODE;

    while ((opt = getopt(argc, argv, MBGC_Params::APPEND_SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                params.forceOverwrite = true;
            break;
            case MBGC_Params::NO_OF_THREADS_OPT:
                newParams.coderThreads = atoi(optarg);
                params.setThreadsLimit(newParams.coderThreads);
                break;
            case MBGC_Params::MATCHER_NO_OF_THREADS_OPT:
                newParams.matcherWorkingThreads = atoi(optarg);
                params.setMatcherThreads(newParams.matcherWorkingThreads);
                break;
            case MBGC_Params::INPUT_FILE_OPT:
                params.setInputFileName(string(optarg));
                break;
            case MBGC_Params::SEQ_CODING_MODE_OPT:
                newParams.setSequentialMatchingMode();
                params.disableLazyDecompression = true;
                break;
#ifdef DEVELOPER_BUILD
            case MBGC_Params::PRINTOUT_PARAMS_OPT:
                printParamsFlag = true;
                break;
            case MBGC_Params::PARAMS_TO_LOGFILE_OPT:
                params.enableLogFile(optarg);
                break;
#endif
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::APP_OUT_TO_CERR_OPT:
                params.setAppOutToCerr();
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printAppendUsage(base_cmd_name, true);
#ifdef DEVELOPER_BUILD
                printAppendDeveloperOptions();
#endif
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n\n");
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
    if (params.inputFileName.empty() && (optind > (argc - 2) || optind < (argc - 3))) {
        fprintf(stderr, "%s: For 'append' command using files list as input expected 2 or 3 arguments after options "
                        "(found %d)\n", base_cmd_name.c_str(), argc - optind);
        cliError = true;
    }
    if (!params.inputFileName.empty() && (optind > (argc - 1) || optind < (argc - 2))) {
        fprintf(stderr, "%s: For 'append' command given single file as input expected 1 or 2 arguments after options "
                        "(found %d)\n", base_cmd_name.c_str(), argc - optind);
        cliError = true;
    }
    if (PgHelpers::numberOfThreads <= 0) {
        cliError = true;
    }
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

    if (params.inputFileName.empty())
        params.setSeqListFileName(argv[optind++]);
    params.setInArchiveFileName(argv[optind++]);
    params.setOutArchiveFileName(optind < argc ? argv[optind++] : params.inArchiveFileName);
    if (params.appOutToCerr) {
        PgHelpers::appout = &std::cerr;
#ifdef DEVELOPER_BUILD
        PgHelpers::devout = &std::cerr;
#endif
    } else if (params.outArchiveFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
        PgHelpers::devout = &null_stream;
    }

    MBGC_Encoder encoder(&params);
    encoder.append(newParams);

    *PgHelpers::appout << "appending time - " << PgHelpers::time_millis() << " [ms]" << endl;

    if (printParamsFlag)
        params.printout();
    return EXIT_SUCCESS;
}


void printCompareUsage(string base_cmd_name, bool details) {
    printVersion(details);
    fprintf(stderr, "\nUsage for counting per-record site differences vs the first (G0) record in archive:\n\t%s"
                    " [-%c <diffBinFile>] <sitesFile> <archiveFile> [<outputFile>]\n", base_cmd_name.c_str(), 'D');
    fprintf(stderr, "Usage for querying a previously dumped diff file (fast, no archive access):\n\t%s"
                    " -%c <diffBinFile> <sitesFile> [<outputFile>]\n", base_cmd_name.c_str(), 'q');
    fprintf(stderr, "Usage for batch query (diff loaded once, N samples):\n\t%s"
                    " -%c -%c <batchListFile> <diffBinFile>\n", base_cmd_name.c_str(), 'q', 'B');
    if (details) {
        fprintf(stderr,
                "\n<sitesFile> name of text file with 1-based site positions (one per line)\n"
                "<archiveFile> mbgc archive filename\n"
                "<diffBinFile> binary file with per-record differing/invalid site lists\n"
                "\t(written with -%c during comparison, queried with -%c)\n"
                "\t(a name ending with .zst is auto zstd-compressed when dumping;\n"
                "\tzstd input is auto-detected and decompressed when querying)\n"
                "<outputFile> output TSV file: record (strain) id, differences count, valid compared sites\n"
                "\t(standard output by default)\n"
                "<batchListFile> text file with one job per line: <sitesFile> <outputFile>\n"
                "\t(blank lines and lines starting with '#' are skipped;\n"
                "\tthe diff file is loaded and parsed only once for all jobs)\n"
                "\nEvery FASTA record in the archive is treated as one strain;\n"
                "the first record of the first file is the reference (G0).\n"
                "Positions where either sequence differs from A/C/G/T (e.g. N or -) are not counted.\n", 'D', 'q');
    }
    fprintf(stderr, "\nBasic options:\n");
    fprintf(stderr, "\t[-%c] overwrite an existing output file\n", MBGC_Params::FORCE_OVERWRITE_OPT);
    printCommonOpts();
}

/* ---------- zstd integration (ported from the mbgc2 wrapper) ----------
 * g -D <diff>.zst : diff is automatically zstd-compressed after comparison
 * g -q <diff>.zst : input is auto-detected (zstd magic) and decompressed before query
 * env: MBGC_ZSTD_LEVEL (default 19), MBGC_ZSTD_THREADS (default 4)          */
static const unsigned char ZSTD_MAGIC[4] = {0x28, 0xB5, 0x2F, 0xFD};

static bool zstd_file_exists(const char *p) { struct stat st; return stat(p, &st) == 0; }

static bool is_zstd_file(const string &p) {
    FILE *f = fopen(p.c_str(), "rb");
    if (!f) return false;
    unsigned char m[4];
    size_t n = fread(m, 1, 4, f);
    fclose(f);
    return n == 4 && memcmp(m, ZSTD_MAGIC, 4) == 0;
}

static bool has_suffix(const string &s, const string &suf) {
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

static int zstd_compress_file(const string &src, const string &dst, int level, int threads) {
    FILE *in = fopen(src.c_str(), "rb");
    FILE *out = fopen(dst.c_str(), "wb");
    if (!in || !out) { if (in) fclose(in); if (out) fclose(out); return -1; }
    ZSTD_CCtx *c = ZSTD_createCCtx();
    if (!c) return -1;
    ZSTD_CCtx_setParameter(c, ZSTD_c_compressionLevel, level);
    ZSTD_CCtx_setParameter(c, ZSTD_c_nbWorkers, threads);
    size_t inCap = 1 << 20;
    size_t outCap = ZSTD_compressBound(inCap);
    void *inBuf = malloc(inCap), *outBuf = malloc(outCap);
    size_t rd;
    while ((rd = fread(inBuf, 1, inCap, in)) > 0) {
        ZSTD_inBuffer ib = {inBuf, rd, 0};
        while (ib.pos < ib.size) {
            ZSTD_outBuffer ob = {outBuf, outCap, 0};
            size_t ret = ZSTD_compressStream2(c, &ob, &ib, ZSTD_e_continue);
            if (ZSTD_isError(ret)) return -1;
            fwrite(outBuf, 1, ob.pos, out);
        }
    }
    ZSTD_inBuffer ib = {NULL, 0, 0};
    size_t ret;
    do {
        ZSTD_outBuffer ob = {outBuf, outCap, 0};
        ret = ZSTD_compressStream2(c, &ob, &ib, ZSTD_e_end);
        if (ZSTD_isError(ret)) return -1;
        fwrite(outBuf, 1, ob.pos, out);
    } while (ret != 0);
    fclose(in); fclose(out);
    ZSTD_freeCCtx(c);
    free(inBuf); free(outBuf);
    return 0;
}

static int zstd_decompress_file(const string &src, const string &dst) {
    FILE *in = fopen(src.c_str(), "rb");
    FILE *out = fopen(dst.c_str(), "wb");
    if (!in || !out) { if (in) fclose(in); if (out) fclose(out); return -1; }
    ZSTD_DCtx *d = ZSTD_createDCtx();
    if (!d) return -1;
    size_t inCap = 1 << 20, outCap = 1 << 22;
    void *inBuf = malloc(inCap), *outBuf = malloc(outCap);
    ZSTD_inBuffer ib = {inBuf, 0, 0};
    int done = 0;
    while (!done) {
        if (ib.pos == ib.size) {
            ib.size = fread(inBuf, 1, inCap, in);
            ib.pos = 0;
        }
        ZSTD_outBuffer ob = {outBuf, outCap, 0};
        size_t ret = ZSTD_decompressStream(d, &ob, &ib);
        if (ZSTD_isError(ret)) return -1;
        fwrite(outBuf, 1, ob.pos, out);
        if (ret == 0) done = 1;
        else if (ib.size == 0 && feof(in)) { return -1; }
    }
    fclose(in); fclose(out);
    ZSTD_freeDCtx(d);
    free(inBuf); free(outBuf);
    return 0;
}

/* read a batch query list: one job per line "<sitesFile> <outputFile>",
 * blank lines and '#' comments are skipped */
static bool loadBatchJobsFile(const string& path, vector<pair<string, string>>& jobs) {
    ifstream in(path);
    if (!in) {
        fprintf(stderr, "ERROR: cannot open batch list file %s\n", path.c_str());
        return false;
    }
    string line;
    uint64_t lineNo = 0;
    auto isWs = [](char c) { return c == ' ' || c == '\t'; };
    while (getline(in, line)) {
        lineNo++;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        size_t i = 0;
        while (i < line.size() && isWs(line[i]))
            i++;
        if (i >= line.size() || line[i] == '#')
            continue;
        size_t j = i;
        while (j < line.size() && !isWs(line[j]))
            j++;
        const string sitesFile = line.substr(i, j - i);
        i = j;
        while (i < line.size() && isWs(line[i]))
            i++;
        size_t k = i;
        while (k < line.size() && !isWs(line[k]))
            k++;
        const string outFile = line.substr(i, k - i);
        size_t t = k;
        while (t < line.size() && isWs(line[t]))
            t++;
        if (sitesFile.empty() || outFile.empty() || t < line.size()) {
            fprintf(stderr, "ERROR: invalid line %llu in batch list %s "
                            "(expected '<sitesFile> <outputFile>')\n",
                    (unsigned long long) lineNo, path.c_str());
            return false;
        }
        jobs.emplace_back(sitesFile, outFile);
    }
    if (jobs.empty()) {
        fprintf(stderr, "ERROR: no jobs found in batch list %s\n", path.c_str());
        return false;
    }
    return true;
}

int compareCommand(std::string base_cmd_name, int argc, char *argv[]) {
    base_cmd_name += string(" ") + string(argv[0]);
    argv[0] = (char *) base_cmd_name.c_str();
    if (argc == 1) {
        printCompareUsage(base_cmd_name, false);
        fprintf(stderr, "\n");
        return EXIT_SUCCESS;
    }
    int opt;
    bool cliError = false;
    bool queryMode = false;
    string dumpFileName;
    string batchListFile;
    MBGC_Params params, inParams;

    while ((opt = getopt(argc, argv, MBGC_Params::COMPARE_SHORT_OPTS.c_str())) != -1) {
        switch (opt) {
            case 'D':
                dumpFileName = optarg;
                break;
            case 'q':
                queryMode = true;
                break;
            case 'B':
                batchListFile = optarg;
                break;
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                params.forceOverwrite = true;
                inParams.forceOverwrite = true;
                break;
            case MBGC_Params::NO_OF_THREADS_OPT:
                params.setThreadsLimit(atoi(optarg));
                inParams.setThreadsLimit(atoi(optarg));
                break;
            case MBGC_Params::IGNORE_FASTA_FILES_PATHS_OPT:
                params.ignoreFastaFilesPath = true;
                inParams.ignoreFastaFilesPath = true;
                break;
            case MBGC_Params::APP_OUT_TO_CERR_OPT:
                params.setAppOutToCerr();
                break;
            case MBGC_Params::VERSION_OPT:
                printVersion(false);
                return EXIT_SUCCESS;
            case MBGC_Params::HELP_OPT:
                printCompareUsage(base_cmd_name, true);
                fprintf(stderr, "\nThe order of all selected options is arbitrary.\n\n");
                return EXIT_SUCCESS;
            default: /* '?' */
                cliError = true;
        }
    }
    if (batchListFile.empty()) {
        if (optind < (argc - 3) || optind > (argc - 2)) {
            fprintf(stderr, "%s: For 'g0-compare' command expected 2 or 3 arguments after options (found %d)\n",
                    base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
    } else {
        if (optind != (argc - 1)) {
            fprintf(stderr, "%s: For 'g0-compare -q -B' command expected 1 argument (diff file) after options (found %d)\n",
                    base_cmd_name.c_str(), argc - optind);
            cliError = true;
        }
    }
    if (!batchListFile.empty() && !queryMode) {
        fprintf(stderr, "%s: -B (batch list) can only be used with -q (query)\n", base_cmd_name.c_str());
        cliError = true;
    }
    if (queryMode && !dumpFileName.empty()) {
        fprintf(stderr, "%s: -D (dump) cannot be combined with -q (query)\n", base_cmd_name.c_str());
        cliError = true;
    }
    if (PgHelpers::numberOfThreads <= 0)
        cliError = true;
    if (cliError) {
        fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
        return EXIT_FAILURE;
    }
    PgHelpers::time_checkpoint();

    if (queryMode) {
        string diffBinFileName = argv[optind++];
        vector<pair<string, string>> jobs;
        if (!batchListFile.empty()) {
            if (!loadBatchJobsFile(batchListFile, jobs))
                return EXIT_FAILURE;
        } else {
            const string sitesFileName = argv[optind++];
            const string outFileName = optind < argc ? argv[optind++] : "";
            if (outFileName.empty() || outFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
                PgHelpers::appout = &null_stream;
#ifdef DEVELOPER_BUILD
                PgHelpers::devout = &null_stream;
#endif
            }
            jobs.emplace_back(sitesFileName, outFileName);
        }
        MBGC_Encoder encoder(&params);
        if (is_zstd_file(diffBinFileName)) {
            /* auto-decompress zstd diff to a temp file once, then query all jobs */
            char tmp[256];
            snprintf(tmp, sizeof(tmp), "/tmp/mbgc_zstd_%d.raw", (int) getpid());
            if (zstd_decompress_file(diffBinFileName, tmp) != 0) {
                fprintf(stderr, "mbgc: zstd decompression failed for %s\n", diffBinFileName.c_str());
                unlink(tmp);
                return EXIT_FAILURE;
            }
            encoder.queryDiffBatch(tmp, jobs, params.forceOverwrite);
            unlink(tmp);
        } else {
            encoder.queryDiffBatch(diffBinFileName, jobs, params.forceOverwrite);
        }
        *PgHelpers::appout << "querying time - " << PgHelpers::time_millis() << " [ms]" << endl;
        return EXIT_SUCCESS;
    }

    const string sitesFileName = argv[optind++];
    inParams.setInArchiveFileName(argv[optind++]);
    const string outFileName = optind < argc ? argv[optind++] : "";
    if (params.appOutToCerr) {
        PgHelpers::appout = &std::cerr;
#ifdef DEVELOPER_BUILD
        PgHelpers::devout = &std::cerr;
#endif
    } else if (outFileName.empty() || outFileName == MBGC_Params::STANDARD_IO_POSIX_ALIAS) {
        PgHelpers::appout = &null_stream;
#ifdef DEVELOPER_BUILD
        PgHelpers::devout = &null_stream;
#endif
    }

    MBGC_Encoder encoder(&params);
    if (!dumpFileName.empty() && has_suffix(dumpFileName, ".zst")) {
        /* zstd-compressed diff: compare writes a temp raw diff, then compress to the target */
        if (zstd_file_exists(dumpFileName.c_str()) && !params.forceOverwrite) {
            fprintf(stderr, "mbgc: %s already exists (use -%c to overwrite)\n",
                    dumpFileName.c_str(), MBGC_Params::FORCE_OVERWRITE_OPT);
            return EXIT_FAILURE;
        }
        char tmp[256];
        snprintf(tmp, sizeof(tmp), "/tmp/mbgc_zstd_%d.raw", (int) getpid());
        encoder.compare(inParams, sitesFileName, outFileName, string(tmp));
        int level = 19, threads = 4;
        const char *e;
        if ((e = getenv("MBGC_ZSTD_LEVEL")) && *e) level = atoi(e);
        if ((e = getenv("MBGC_ZSTD_THREADS")) && *e) threads = atoi(e);
        if (zstd_compress_file(tmp, dumpFileName, level, threads) != 0) {
            fprintf(stderr, "mbgc: zstd compression failed for %s\n", dumpFileName.c_str());
            unlink(tmp);
            return EXIT_FAILURE;
        }
        unlink(tmp);
        *PgHelpers::appout << "diff compressed (zstd level " << level << ") -> " << dumpFileName << endl;
    } else {
        encoder.compare(inParams, sitesFileName, outFileName, dumpFileName);
    }

    *PgHelpers::appout << "comparing time - " << PgHelpers::time_millis() << " [ms]" << endl;
    return EXIT_SUCCESS;
}

void printPreStrainUsage(string base_cmd_name, bool details) {
    printVersion(details);
    fprintf(stderr, "\nUsage for preStrain representative pipeline (sample scoring + selection + extraction, one step):\n");
    fprintf(stderr, "\t%s [-%c] [-%c <ref_num>] [-d <DBdir>] <diff.zst|diff.raw> <cluster.csv> <archive.mbgc> <pickX> <outDir>\n",
            base_cmd_name.c_str(), MBGC_Params::FORCE_OVERWRITE_OPT, 'n');
    fprintf(stderr, "\t  or: %s [-%c] [-%c <ref_num>] [-d <DBdir>] <really-ref.fa> <ref-coords.txt> <diff.zst|diff.raw> <cluster.csv> <archive.mbgc> <pickX> <outDir>\n",
            base_cmd_name.c_str(), MBGC_Params::FORCE_OVERWRITE_OPT, 'n');
    fprintf(stderr, "\t  or: %s [-%c] [-%c <ref_num>] <G0.fa> <really-ref.fa> <ref-coords.txt> <diff.zst|diff.raw> <cluster.csv> <archive.mbgc> <pickX> <outDir>\n",
            base_cmd_name.c_str(), MBGC_Params::FORCE_OVERWRITE_OPT, 'n');
    fprintf(stderr, "\t  BATCH: %s -%c <DBdir> <samples.list> <pickX>\n",
            base_cmd_name.c_str(), 'd');
    fprintf(stderr, "\nBatch mode: <samples.list> has one SampleName per line (first column; the\n");
    fprintf(stderr, "rest is ignored). For every DB listed in <DBdir>/DB.list, the diff is loaded\n");
    fprintf(stderr, "exactly once and all samples are scored in a single pass. Each sample reads\n");
    fprintf(stderr, "YC-result/<SampleName>/output/{really-ref.fa,ref-coords.txt}; results go to\n");
    fprintf(stderr, "YC-result/<SampleName>/output/rep_<db>.tsv and .../rep_<db>/.\n");
    fprintf(stderr, "\nIn the 5-argument form, <really-ref.fa> and <ref-coords.txt> default to\n");
    fprintf(stderr, "./output/really-ref.fa and ./output/ref-coords.txt (current directory).\n");
    fprintf(stderr, "-%c <ref_num>: pre-computed sample score vs G0; skips step 1. The score depends\n",
            'n');
    fprintf(stderr, "only on the sample and G0, not on the DB - compute it once per sample and\n");
    fprintf(stderr, "reuse it for every DB with -%c.\n", 'n');
    fprintf(stderr, "G0 resolution order: explicit <G0.fa> argument, then -d <DBdir>/G0.fa,\n");
    fprintf(stderr, "then the archive's first strain (as recorded in the diff file header).\n");
    fprintf(stderr, "<really-ref.fa>  sample FASTA to score against G0 at the coords sites\n");
    fprintf(stderr, "<ref-coords.txt> 1-based site positions, one per line (same file used for the DB score)\n");
    fprintf(stderr, "<diff.zst>       compressed diff DB (output of 'g -D')\n");
    fprintf(stderr, "<cluster.csv>    strain,cluster_id mapping (one per line)\n");
    fprintf(stderr, "<archive.mbgc>   compressed archive from which representative strains are extracted\n");
    fprintf(stderr, "<pickX>          number of representative strains to pick\n");
    fprintf(stderr, "<outDir>         output directory for the extracted representative FASTA files\n");
    if (details) {
        fprintf(stderr, "\nPipeline (equivalent of the 4 original steps, run here in memory):\n");
        fprintf(stderr, "\t1. count sample diffs vs G0 at coords sites (g0_diff.pl ref_num)\n");
        fprintf(stderr, "\t2. score every DB strain at the same sites (g -q)\n");
        fprintf(stderr, "\t3. per cluster pick the closest strain, then top <pickX> globally (cal2.sh)\n");
        fprintf(stderr, "\t4. extract the picked strains from the archive (d -E)\n");
    }
}

// read the G0 record name from the mDIF diff header (little-endian)
static string readDiffG0Name(const string& rawDiffPath) {
    ifstream in(rawDiffPath, ios::binary);
    if (!in) {
        fprintf(stderr, "ERROR: cannot open %s\n", rawDiffPath.c_str());
        exit(EXIT_FAILURE);
    }
    char magic[4];
    in.read(magic, 4);
    if (memcmp(magic, "mDIF", 4) != 0) {
        fprintf(stderr, "ERROR: %s is not a valid diff file (bad magic).\n", rawDiffPath.c_str());
        exit(EXIT_FAILURE);
    }
    uint32_t version;
    uint64_t g0len;
    uint16_t nameLen;
    in.read((char*) &version, 4);
    in.read((char*) &g0len, 8);
    in.read((char*) &nameLen, 2);
    string name(nameLen, '\0');
    if (nameLen)
        in.read(&name[0], nameLen);
    return name;
}

// ---- per-sample extraction helpers: really-ref scoring + poc-site slicing ----

// read a single-record FASTA, return the uppercased sequence (no '>' header)
static string readFastaSeq(const string& path) {
    ifstream in(path);
    if (!in) return string();
    string line, seq;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '>') continue;
        for (char c : line) {
            if (c == '\r') continue;
            seq.push_back((char) toupper((unsigned char) c));
        }
    }
    return seq;
}

// load 1-based positions (poc.txt), one integer per line
static vector<uint32_t> loadPoc(const string& path) {
    vector<uint32_t> poc;
    ifstream in(path);
    uint32_t x;
    while (in >> x) poc.push_back(x);
    return poc;
}

// SNP distance: count positions where both bases are called (A/C/G/T) and differ;
// N/. /- gaps are ignored on either side
static uint32_t snpDist(const string& a, const string& b) {
    size_t n = min(a.size(), b.size());
    uint32_t d = 0;
    for (size_t i = 0; i < n; i++) {
        char x = a[i], y = b[i];
        if (x != y && x != 'N' && x != '.' && x != '-' && y != 'N' && y != '.' && y != '-')
            d++;
    }
    return d;
}

// keep only the bases at the poc sites (1-based positions of really-ref coords)
static string sliceByPoc(const string& seq, const vector<uint32_t>& poc) {
    string out;
    out.reserve(poc.size());
    for (uint32_t p : poc) {
        size_t idx = (size_t) p - 1; // 1-based -> 0-based
        out.push_back(idx < seq.size() ? seq[idx] : 'N');
    }
    return out;
}

// recursively remove a directory tree (files + subdirs)
static void removeDirRecursive(const string& path) {
    DIR* d = opendir(path.c_str());
    if (!d) return;
    struct dirent* e;
    while ((e = readdir(d)) != nullptr) {
        if (strcmp(e->d_name, ".") == 0 || strcmp(e->d_name, "..") == 0) continue;
        string p = path + "/" + e->d_name;
        struct stat st;
        if (lstat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode))
            removeDirRecursive(p);
        else
            unlink(p.c_str());
    }
    closedir(d);
    rmdir(path.c_str());
}

// extract a single strain (exact record name, d -E semantics) from the archive
// into outDir and return the path of the extracted FASTA file
static string extractStrainFasta(const string& archiveFile, const string& strainName,
                                 const string& outDir, bool forceOverwrite) {
    string patFile = "/tmp/mbgc_ps_g0pat_" + to_string((long) getpid()) + ".txt";
    {
        ofstream pf(patFile);
        pf << strainName << "\n";
    }
    MBGC_Params dParams;
    dParams.forceOverwrite = forceOverwrite;
    dParams.setInArchiveFileName(archiveFile);
    dParams.setOutputPath(outDir);
    dParams.setInputFileName(patFile);
    ostream* savedAppout = PgHelpers::appout;
    PgHelpers::appout = &cerr;
    MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&dParams, string());
    decoder->decode();
    delete(decoder);
    PgHelpers::appout = savedAppout;
    unlink(patFile.c_str());

    DIR* dir = opendir(outDir.c_str());
    if (!dir) {
        fprintf(stderr, "ERROR: cannot open extraction dir %s\n", outDir.c_str());
        exit(EXIT_FAILURE);
    }
    struct dirent* ent;
    string fastaPath;
    while ((ent = readdir(dir)) != nullptr) {
        if (ent->d_name[0] == '.')
            continue;
        fastaPath = outDir + "/" + ent->d_name;
        break;
    }
    closedir(dir);
    if (fastaPath.empty()) {
        fprintf(stderr, "ERROR: strain '%s' not found in archive %s\n", strainName.c_str(), archiveFile.c_str());
        exit(EXIT_FAILURE);
    }
    return fastaPath;
}

int preStrainCommand(string base_cmd_name, int argc, char *argv[]) {
    base_cmd_name += string(" ") + string(argv[0]);
    argv[0] = (char *) base_cmd_name.c_str();
    if (argc == 1) {
        printPreStrainUsage(base_cmd_name, false);
        fprintf(stderr, "\n");
        return EXIT_SUCCESS;
    }
    bool forceOverwrite = false;
    string dbDir;
    uint64_t refNumArg = 0;
    bool refNumProvided = false;
    MBGC_Params params;
    int opt;
    while ((opt = getopt(argc, argv, "fd:n:")) != -1) {
        switch (opt) {
            case MBGC_Params::FORCE_OVERWRITE_OPT:
                forceOverwrite = true;
                break;
            case 'd': // database directory: G0.fa lives at <dbDir>/G0.fa
                dbDir = optarg;
                break;
            case 'n': // pre-computed sample score vs G0 (skips step 1; sample-specific, DB-independent)
                refNumArg = strtoull(optarg, nullptr, 10);
                refNumProvided = true;
                break;
            default:
                fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
                return EXIT_FAILURE;
        }
    }
    // four forms:
    //   3 args:  <samples.list> <pickX> [max_dist]   (batch mode; requires -d <DBdir>; all DBs in DB.list)
    //   2 args:  <samples.list> <pickX>   (batch mode without distance filter)
    //   8 args: <G0.fa> <ref.fa> <coords.txt> <diff> <clu> <mbgc> <pickX> <outDir>
    //   7 args: <ref.fa> <coords.txt> <diff> <clu> <mbgc> <pickX> <outDir>
    //   5 args: <diff> <clu> <mbgc> <pickX> <outDir>  (ref/coords default to ./output/)
    const int nargs = argc - optind;
    if (nargs != 8 && nargs != 7 && nargs != 5 && nargs != 2 && nargs != 3) {
        fprintf(stderr, "%s: For 'preStrain' command expected 2, 3, 5, 7 or 8 arguments after options (found %d)\n",
                base_cmd_name.c_str(), nargs);
        fprintf(stderr, "try '%s -%c' for more information\n", base_cmd_name.c_str(), MBGC_Params::HELP_OPT);
        return EXIT_FAILURE;
    }

    // ---- batch mode: <samples.list> <pickX> [max_dist] with -d <DBdir> ----
    if (nargs == 2 || nargs == 3) {
        if (dbDir.empty()) {
            fprintf(stderr, "%s: batch mode requires -d <DBdir>\n", base_cmd_name.c_str());
            return EXIT_FAILURE;
        }
        const string samplesListFile = argv[optind];
        const uint64_t pickXB = strtoull(argv[optind + 1], nullptr, 10);
        if (pickXB == 0) {
            fprintf(stderr, "%s: <pickX> must be a positive integer\n", base_cmd_name.c_str());
            return EXIT_FAILURE;
        }
        // optional max_dist: strains with really-ref distance > max_dist are skipped
        // (no FASTA written, no tsv row); 0 or absent = no filter
        const uint64_t maxDistB = (nargs == 3) ? strtoull(argv[optind + 2], nullptr, 10) : 0;
        if (maxDistB > 0)
            fprintf(stderr, "batch distance filter: keep strains with really-ref dist <= %llu\n",
                    (unsigned long long) maxDistB);
        // sample names: first column of each non-comment line
        vector<string> sampleNames;
        {
            ifstream sl(samplesListFile);
            if (!sl) {
                fprintf(stderr, "%s: cannot open samples list %s\n", base_cmd_name.c_str(), samplesListFile.c_str());
                return EXIT_FAILURE;
            }
            string line;
            while (getline(sl, line)) {
                if (!line.empty() && line.back() == '\r')
                    line.pop_back();
                size_t b = line.find_first_not_of(" \t");
                if (b == string::npos || line[b] == '#')
                    continue;
                size_t e = line.find_first_of(" \t", b);
                sampleNames.push_back(line.substr(b, e == string::npos ? string::npos : e - b));
            }
        }
        if (sampleNames.empty()) {
            fprintf(stderr, "%s: no samples found in %s\n", base_cmd_name.c_str(), samplesListFile.c_str());
            return EXIT_FAILURE;
        }
        // DB names: from <DBdir>/DB.list if present, else scan DB/*.zst
        vector<string> dbs;
        {
            const string listPath = dbDir + "/DB.list";
            ifstream dbf(listPath);
            if (dbf) {
                string line;
                while (getline(dbf, line)) {
                    if (!line.empty() && line.back() == '\r')
                        line.pop_back();
                    size_t b = line.find_first_not_of(" \t");
                    if (b == string::npos || line[b] == '#')
                        continue;
                    size_t e = line.find_last_not_of(" \t");
                    dbs.push_back(line.substr(b, e - b + 1));
                }
            } else {
                // auto-discover: list *.zst in DB/
                const string dbSub = dbDir + "/DB";
                DIR* dir = opendir(dbSub.c_str());
                if (!dir) {
                    fprintf(stderr, "%s: cannot open %s/DB (no DB.list either)\n", base_cmd_name.c_str(), dbDir.c_str());
                    return EXIT_FAILURE;
                }
                struct dirent* ent;
                while ((ent = readdir(dir)) != nullptr) {
                    string name = ent->d_name;
                    if (name.size() > 4 && name.compare(name.size() - 4, 4, ".zst") == 0)
                        dbs.push_back(name.substr(0, name.size() - 4));
                }
                closedir(dir);
                sort(dbs.begin(), dbs.end());
            }
        }
        const string g0File = dbDir + "/G0.fa";
        {
            ifstream gf(g0File);
            if (!gf) {
                fprintf(stderr, "%s: G0 not found: %s\n", base_cmd_name.c_str(), g0File.c_str());
                return EXIT_FAILURE;
            }
        }
        fprintf(stderr, "batch: %zu samples, %zu DBs, G0=%s\n", sampleNames.size(), dbs.size(), g0File.c_str());

        MBGC_Encoder encoder(&params);
        for (const string& db : dbs) {
            const string zstF = dbDir + "/DB/" + db + ".zst";
            const string cluF = dbDir + "/DB/" + db + ".clu";
            const string mbgcF = dbDir + "/DB/" + db + ".mbgc";
            if (access(zstF.c_str(), R_OK) != 0 || access(cluF.c_str(), R_OK) != 0 || access(mbgcF.c_str(), R_OK) != 0) {
                fprintf(stderr, "SKIP db %s: incomplete files\n", db.c_str());
                continue;
            }
            // samples with ready output/ files
            vector<MBGC_Encoder::PreSample> jobs;
            for (const string& name : sampleNames) {
                const string refF = "YC-result/" + name + "/output/really-ref.fa";
                const string coF = "YC-result/" + name + "/output/ref-coords.txt";
                if (access(refF.c_str(), R_OK) != 0 || access(coF.c_str(), R_OK) != 0) {
                    fprintf(stderr, "SKIP sample %s: output files missing\n", name.c_str());
                    continue;
                }
                MBGC_Encoder::PreSample ps;
                ps.name = name;
                ps.refFile = refF;
                ps.coordsFile = coF;
                jobs.push_back(ps);
            }
            if (jobs.empty()) {
                fprintf(stderr, "no valid samples for db %s\n", db.c_str());
                continue;
            }
            // decompress the diff .zst once for all samples of this db
            string diffRaw = zstF;
            string tmpRaw;
            if (is_zstd_file(zstF)) {
                char tmp[256];
                snprintf(tmp, sizeof(tmp), "/tmp/mbgc_psb_%d_%s.raw", (int) getpid(), db.c_str());
                if (zstd_decompress_file(zstF, tmp) != 0) {
                    fprintf(stderr, "zstd decompression failed for %s\n", zstF.c_str());
                    continue;
                }
                diffRaw = tmp;
                tmpRaw = tmp;
            }
            fprintf(stderr, "=== DB %s: scoring %zu samples with a single diff load ===\n", db.c_str(), jobs.size());
            auto tScore0 = std::chrono::steady_clock::now();
            vector<MBGC_Encoder::PreSampleResult> res = encoder.preStrainBatch(g0File, jobs, diffRaw, cluF, pickXB);
            auto tScore1 = std::chrono::steady_clock::now();
            if (!tmpRaw.empty())
                unlink(tmpRaw.c_str());

            // ---- per-sample extraction: decode once to a temp dir, then for each
            //      sample score every rep vs its own really-ref.fa and emit a sliced
            //      FASTA at the sample's poc.txt sites, writing to the sample's own
            //      output dir (YC-result/<sample>/output/fa_<db>/) ----
            // temp dir on /tmp (native disk), not /mnt/c 9P
            const string tmpDir = "/tmp/mbgc_psb_" + to_string((long) getpid()) + "_" + db;
            removeDirRecursive(tmpDir);
            mkdir(tmpDir.c_str(), 0755);

            std::set<string> unionNames;
            for (const auto& r : res)
                for (const string& n : r.reps)
                    unionNames.insert(n);

            const string namesFile = "/tmp/mbgc_psb_names_" + to_string((long) getpid()) + "_" + db + ".txt";
            {
                ofstream nf(namesFile);
                for (const string& n : unionNames)
                    nf << n << "\n";
            }
            auto tDec0 = std::chrono::steady_clock::now();
            {
                MBGC_Params dParams;
                dParams.forceOverwrite = true;
                dParams.setInArchiveFileName(mbgcF);
                dParams.setOutputPath(tmpDir);
                dParams.setInputFileName(namesFile);
                ostream* savedAppout = PgHelpers::appout;
                PgHelpers::appout = &cerr;
                MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&dParams, string());
                decoder->decode();
                delete(decoder);
                PgHelpers::appout = savedAppout;
            }
            unlink(namesFile.c_str());
            auto tDec1 = std::chrono::steady_clock::now();

            // load per-sample context: really-ref sequence + poc sites
            vector<pair<string, vector<uint32_t>>> ctx(res.size());
            for (size_t si = 0; si < res.size(); si++) {
                ctx[si].first = readFastaSeq("YC-result/" + res[si].name + "/output/really-ref.fa");
                ctx[si].second = loadPoc("YC-result/" + res[si].name + "/output/poc.txt");
                if (ctx[si].first.empty())
                    fprintf(stderr, "WARN: really-ref empty for %s\n", res[si].name.c_str());
            }

            auto tPost0 = std::chrono::steady_clock::now();
            // samples are independent: process in parallel (each writes its own dir)
            omp_set_num_threads(4);
            #pragma omp parallel for schedule(dynamic)
            for (int si = 0; si < (int) res.size(); si++) {
                const string outBase = "YC-result/" + res[si].name + "/output/";
                const string faDir = outBase + "fa_" + db;
                removeDirRecursive(faDir);
                mkdir(faDir.c_str(), 0755);

                // keep the original cluster id per strain from the selection repText
                map<string, string> cidMap;
                {
                    istringstream iss(res[si].repText);
                    string st, cid;
                    unsigned long long av;
                    while (iss >> st >> av >> cid)
                        cidMap[st] = cid;
                }

                string repTsv;
                size_t nKept = 0, nFiltered = 0;
                for (const string& n : res[si].reps) {
                    string seq = readFastaSeq(tmpDir + "/" + n);
                    if (seq.empty()) {
                        fprintf(stderr, "WARN: no sequence extracted for %s\n", n.c_str());
                        continue;
                    }
                    uint32_t d = snpDist(seq, ctx[si].first);
                    // distance filter: skip strains too far from this sample's really-ref
                    if (maxDistB > 0 && d > maxDistB) {
                        nFiltered++;
                        continue;
                    }
                    string sliced = sliceByPoc(seq, ctx[si].second);
                    {
                        ofstream of(faDir + "/" + n);
                        of << ">" << n << "\n";
                        for (size_t i = 0; i < sliced.size(); i += 80)
                            of << sliced.substr(i, 80) << "\n";
                    }
                    repTsv += n + "\t" + to_string(d) + "\t" + cidMap[n] + "\n";
                    nKept++;
                }
                {
                    ofstream tf(outBase + "rep_" + db + ".tsv");
                    tf << repTsv;
                }
                fprintf(stderr, "sample %s / db %s: ref_num=%llu, %zu reps scored vs really-ref, kept=%zu filtered(dist>%llu)=%zu, %zu poc sites -> %s\n",
                        res[si].name.c_str(), db.c_str(), (unsigned long long) res[si].refNum,
                        res[si].reps.size(), nKept, (unsigned long long) maxDistB, nFiltered,
                        ctx[si].second.size(), faDir.c_str());
            }
            removeDirRecursive(tmpDir);
            auto tPost1 = std::chrono::steady_clock::now();
            fprintf(stderr, "DB %s: scoring %.2fs, decode %.2fs, per-sample score+sliced-extract %.2fs, %zu reps x %zu samples\n",
                    db.c_str(),
                    std::chrono::duration<double>(tScore1 - tScore0).count(),
                    std::chrono::duration<double>(tDec1 - tDec0).count(),
                    std::chrono::duration<double>(tPost1 - tPost0).count(),
                    unionNames.size(), res.size());
        }
        return EXIT_SUCCESS;
    }

    string g0File, refFile, coordsFile;
    int o = optind;
    if (nargs == 8) {
        g0File = argv[o++];
        refFile = argv[o++];
        coordsFile = argv[o++];
    } else if (nargs == 7) {
        refFile = argv[o++];
        coordsFile = argv[o++];
    }
    const string diffFile = argv[o];
    const string clusterFile = argv[o + 1];
    const string archiveFile = argv[o + 2];
    const uint64_t pickX = strtoull(argv[o + 3], nullptr, 10);
    const string outDir = argv[o + 4];
    if (nargs == 5) {
        refFile = "./output/really-ref.fa";
        coordsFile = "./output/ref-coords.txt";
    }
    if (pickX == 0) {
        fprintf(stderr, "%s: <pickX> must be a positive integer\n", base_cmd_name.c_str());
        return EXIT_FAILURE;
    }
    {
        ifstream rf(refFile);
        if (!rf) {
            fprintf(stderr, "%s: cannot open reference FASTA %s\n", base_cmd_name.c_str(), refFile.c_str());
            return EXIT_FAILURE;
        }
        ifstream cf(coordsFile);
        if (!cf) {
            fprintf(stderr, "%s: cannot open coords file %s\n", base_cmd_name.c_str(), coordsFile.c_str());
            return EXIT_FAILURE;
        }
    }

    // diff DB: transparently decompress a .zst diff to a temp raw file (loaded once)
    string diffRaw = diffFile;
    string tmpRaw;
    if (is_zstd_file(diffFile)) {
        char tmp[256];
        snprintf(tmp, sizeof(tmp), "/tmp/mbgc_ps_%d.raw", (int) getpid());
        if (zstd_decompress_file(diffFile, tmp) != 0) {
            fprintf(stderr, "mbgc: zstd decompression failed for %s\n", diffFile.c_str());
            unlink(tmp);
            return EXIT_FAILURE;
        }
        diffRaw = tmp;
        tmpRaw = tmp;
    }

    // G0 resolution order: explicit <G0.fa> arg > -d <DBdir>/G0.fa > first strain of the archive
    string g0TmpFile;
    if (g0File.empty() && !dbDir.empty()) {
        const string g0FromDir = dbDir + "/G0.fa";
        if (access(g0FromDir.c_str(), R_OK) == 0) {
            g0File = g0FromDir;
            fprintf(stderr, "G0 (from -d dir): %s\n", g0File.c_str());
        }
    }
    if (g0File.empty()) {
        const string g0name = readDiffG0Name(diffRaw);
        const string g0Dir = "/tmp/mbgc_ps_g0_" + to_string((long) getpid());
        mkdir(g0Dir.c_str(), 0755);
        g0File = extractStrainFasta(archiveFile, g0name, g0Dir, true);
        g0TmpFile = g0File;
        fprintf(stderr, "G0 (first strain of archive): %s\n", g0name.c_str());
    }

    MBGC_Encoder encoder(&params);
    string repText;
    vector<string> names;
    names = encoder.preStrainSelect(g0File, refFile, coordsFile, diffRaw, clusterFile, pickX,
                                    refNumArg, refNumProvided, repText);
    if (!tmpRaw.empty())
        unlink(tmpRaw.c_str());
    if (!g0TmpFile.empty())
        unlink(g0TmpFile.c_str());
    if (names.empty()) {
        fprintf(stderr, "preStrain: no representative strains selected\n");
        return EXIT_FAILURE;
    }

    // extraction: reuse the archive decoder with a patterns list (d -E semantics)
    string namesFile = "/tmp/mbgc_ps_names_" + to_string((long) getpid()) + ".txt";
    {
        ofstream nf(namesFile);
        for (const string& n : names)
            nf << n << "\n";
    }
    MBGC_Params dParams;
    dParams.forceOverwrite = forceOverwrite;
    dParams.setInArchiveFileName(archiveFile);
    dParams.setOutputPath(outDir);
    dParams.setInputFileName(namesFile);
    // keep the rep table clean: decoder progress messages go to stderr here
    ostream* savedAppout = PgHelpers::appout;
    PgHelpers::appout = &cerr;
    MGMP_Decoder_API* decoder = MBGC_Decoder<true>::getInstance(&dParams, string());
    decoder->decode();
    delete(decoder);
    PgHelpers::appout = savedAppout;
    unlink(namesFile.c_str());

    *PgHelpers::appout << repText;
    fprintf(stderr, "preStrain: selected %zu representatives -> %s\n", names.size(), outDir.c_str());
    return EXIT_SUCCESS;
}

void printCommands(string base_toolname) {
    fprintf(stderr, "\nUsage:\t%s <command> [options]\n", base_toolname.c_str());
    fprintf(stderr, "\nAvailable commands (%c - default):\n", MBGC_Params::INFO_CMD);
    for(int i = 0; i < MBGC_Params::cliCommandsCount; i++) {
        fprintf(stderr, "\t%c\t%s\n",
                MBGC_Params::CLI_COMMANDS[i][0], MBGC_Params::CLI_COMMANDS[i]);
    }
}

int printCommandsUsage(string base_toolname, string unknown_command_arg = "") {
    bool isCommandWrong = !unknown_command_arg.empty();
    if (!isCommandWrong)
        printVersion(false);
    else
        fprintf(stderr, "%s: invalid command -- '%s'\n", base_toolname.c_str(), unknown_command_arg.c_str());
    printCommands(base_toolname);
    if (!isCommandWrong) {
        fprintf(stderr, "\nCommon options:\n");
        printCommonOpts();
#ifdef DEVELOPER_BUILD
        fprintf(stderr, "\n------------------ DEVELOPER OPTIONS ----------------\n");
        printCommonDevOpts();
#endif
    }
    fprintf(stderr, "\n");
    return isCommandWrong ? EXIT_FAILURE : EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
    std::string path = argv[0];
    std::string base_cmd_name = path.substr(path.find_last_of("/\\") + 1);
    if (argc == 1)
        return printCommandsUsage(base_cmd_name);
#ifndef DEVELOPER_BUILD
    PgHelpers::devout = &null_stream;
#endif
    string command_arg = string(argv[1]);
    if (command_arg.size() > 1)
        return decodeCommands(base_cmd_name, argc, argv);

    switch (command_arg[0]) {
        case MBGC_Params::COMPRESS_CMD:
        case MBGC_Params::REPACK_CMD:
            return encodeCommand(base_cmd_name, argc - 1, argv + 1);
        case MBGC_Params::DECOMPRESS_CMD:
        case MBGC_Params::INFO_CMD:
            return decodeCommands(base_cmd_name, argc - 1, argv + 1);
        case MBGC_Params::APPEND_CMD:
            return appendCommand(base_cmd_name, argc - 1, argv + 1);
        case MBGC_Params::COMPARE_CMD:
            return compareCommand(base_cmd_name, argc - 1, argv + 1);
        case MBGC_Params::PRESTRAIN_CMD:
            return preStrainCommand(base_cmd_name, argc - 1, argv + 1);
#ifdef DEVELOPER_BUILD
        case MBGC_Params::VALIDATE_CMD:
            return decodeCommands(base_cmd_name, argc - 1, argv + 1);
#endif
        default: /* '?' */
            printCommandsUsage(base_cmd_name, command_arg);
            return EXIT_FAILURE;
    }
}

