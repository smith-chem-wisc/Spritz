configfile: "config/config.yaml"
report: "report/workflow.rst"
include: "rules/common.smk"

rule all:
    input: all_output # see common.smk

rule setup:
    input: setup_output # see common.smk
    output: "../resources/setup.txt"
    log: "../resources/setup.log"
    conda: "envs/default.yaml"
    shell: "touch {output}"

rule clean:
    log: expand("{dir}/clean.log", dir=config['analysis_directory'])
    conda: "envs/proteogenomics.yaml"
    shell:
        "rm -rf ../resources/ ptmlist.txt PSI-MOD.obo.xml &&"
        "cd ../TransferUniProtModifications && dotnet clean && cd .. 2> {log}"

rule prose:
    output: "{dir}/prose.txt"
    log: "{dir}/prose.log"
    conda: "envs/default.yaml"
    shell: "python scripts/prose.py {output} 2> {log}"

include: "rules/downloads.smk"
include: "rules/align.smk"
include: "rules/variants.smk"
include: "rules/isoforms.smk"
include: "rules/proteogenomics.smk"
include: "rules/quant.smk"
