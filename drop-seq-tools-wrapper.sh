#!/bin/bash
#Author: Paul Hoffman
set -eo pipefail

JAR_DIR="$(dirname $0)/jar"

declare -a PROGRAMS=(
    'BAMTagHistogram'
    'BAMTagofTagCounts'
    'BaseDistributionAtReadPosition'
    'CollapseBarcodesInPlace'
    'CollapseTagWithContext'
    'ConvertToRefFlat'
    'CreateIntervalsFiles'
    'DetectBeadSynthesisErrors'
    'DigitalExpression'
    'FilterBAM'
    'FilterBAMByTag'
    'GatherGeneGCLength'
    'GatherMolecularBarcodeDistributionByGene'
    'GatherReadQualityMetrics'
    'PolyATrimmer'
    'ReduceGTF'
    'SelectCellsByNumTranscripts'
    'SingleCellRnaSeqMetricsCollector'
    'TagBamWithReadSequenceExtended'
    'TagReadWithGeneExon'
    'TagReadWithInterval'
    'TrimStartingSequence'
    'ValidateReference'
)

function Usage() {
    echo -e"\
Usage: $(basname $0) -p|--program <program> [-m|--memory <memory>] [-t|--temp-dir <temp>] [program_args [program_args...]]
" >&2
    exit 1
}

export -f Usage

declare -a ARGUMENTS=()

[[ $# -lt 1 ]] && Usage

while [[ $# -ge 1 ]]; do
    KEY=$1
    case $KEY in
        -p|--program)
            PROGRAM=$2
            shift
            ;;
        -m|--memory)
            XMX=$2
            shift
            ;;
        -t|--temp-dir)
            TEMP=$2
            shift
            ;;
        *)
            ARGUMENTS+=($1)
            ;;
    esac
    shift
done

[[ -z "${XMX}" ]] && MEM='2g'
[[ -z "${TEMP}" ]] && TEMP='/tmp'
[[ -z "${PROGRAM}" ]] && (echo "'-p|--program' is missing" >&2; exit 1)

[[ ${PROGRAMS[@]} =~ ${PROGRAM} ]] || (echo "Invalid program: ${PROGRAM}" >&2; exit 1)

[[ -d "${TEMP}" ]] || (set -x; mkdir -p "${TEMP}")
(set -x; java -Xmx${XMX} -Djava.io.tmpdir=${TEMP} -jar "${JAR_DIR}/dropseq.jar" "${PROGRAM}" ${ARGUMENTS[@]})
