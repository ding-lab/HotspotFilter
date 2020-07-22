class: CommandLineTool
cwlVersion: v1.0
id: hotspotfilter
baseCommand:
  - /bin/bash
  - /opt/HotspotFilter/src/hotspot_filter.sh
inputs:
  - id: VCF_A
    type: File
    inputBinding:
      position: 0
      prefix: '-A'
  - id: VCF_B
    type: File?
    inputBinding:
      position: 0
      prefix: '-B'
  - id: BED
    type: File
    inputBinding:
      position: 0
      prefix: '-D'
  - id: keep_only_pass
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-R'
    doc: Retain only variants with FILTER value of PASS or .
  - id: retain_all
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-r'
    doc: Retain all variants in A, mark exclusions with VCF FILTER field. Ignored if VCF_B specified
  - id: filter_name
    type: string?
    inputBinding:
      position: 0
      prefix: '-F'
    doc: specify filter name in VCF file.  Default is 'hotspot'
outputs:
  - id: output
    type: File
    outputBinding:
      glob: output/HotspotFiltered.vcf
label: HotspotFilter
arguments:
  - position: 0
    prefix: '-o'
    valueFrom: output/HotspotFiltered.vcf
requirements:
  - class: ResourceRequirement
    ramMin: 1000
  - class: DockerRequirement
    dockerPull: 'mwyczalkowski/hotspot_filter:20200617'
