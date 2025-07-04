#!/bin/bash -euo pipefail
workflow-glue check_bam_headers_in_dir input_dir > env.vars
source env.vars

# capture process environment
set +u
set +e
cd "$NXF_TASK_WORKDIR"

nxf_eval_cmd() {
    {
        IFS=$'\n' read -r -d '' "${1}";
        IFS=$'\n' read -r -d '' "${2}";
        (IFS=$'\n' read -r -d '' _ERRNO_; return ${_ERRNO_});
    } < <((printf '\0%s\0%d\0' "$(((({ shift 2; "${@}"; echo "${?}" 1>&3-; } | tr -d '\0' 1>&4-) 4>&2- 2>&1- | tr -d '\0' 1>&4-) 3>&1- | exit "$(cat)") 4>&1-)" "${?}" 1>&2) 2>&1)
}

echo '' > .command.env
#
echo IS_UNALIGNED="${IS_UNALIGNED[@]}" >> .command.env
echo /IS_UNALIGNED/ >> .command.env
#
echo MIXED_HEADERS="${MIXED_HEADERS[@]}" >> .command.env
echo /MIXED_HEADERS/ >> .command.env
#
echo IS_SORTED="${IS_SORTED[@]}" >> .command.env
echo /IS_SORTED/ >> .command.env
