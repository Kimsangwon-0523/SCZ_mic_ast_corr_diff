# MCM perturbation 분석

결과 노트북은 GSE282665와 GSE228927 astrocyte의 ACHES 및 SNAP-a를 보여준다. 두 코호트는 [그림 24](../../figures/24.pdf)와 [그림 25](../../figures/25.pdf)에 각각 제시한다.

- [00_inspect_perturbation_data.ipynb](00_inspect_perturbation_data.ipynb)에서 발현 파일과 원본 열을 확인한다.
- [01_review_program_scores.ipynb](01_review_program_scores.ipynb)에서 데이터셋별 표본 점수, 비교 괄호의 별표와 q값, 조건 간 차이의 95% 신뢰구간을 확인한다.

## 다시 계산하기

GEO supplementary files에서 발현 자료를 내려받아 아래 경로에 둔다. 압축된 파일은 푼 뒤 원래 파일명을 유지한다. 이 발현 파일들은 Git 추적에서 제외한다.

| 출처 | 로컬 폴더 | 점수 계산 입력 |
| --- | --- | --- |
| [GSE228927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228927) | `data/perturbation_validation/GSE228927/` | `GSE228927_FPKM_anno.xls` (탭 구분 텍스트) |
| [GSE282665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282665) | `data/perturbation_validation/GSE282665/` | `GSE282665_astrocyte.GSE282665_microglia.media.salmon.merged.gene_tpm.tsv` 및 `GSE282665_motor.neuron.GSE282665_microglia.media.salmon.merged.gene_tpm.tsv` |

탐색 노트북은 GSE282665의 `GSE282665_microglia.cxcr4.gpnmb.sirna.salmon.merged.gene_tpm.tsv`도 읽는다.

저장소 루트에서 pandas, numpy, scipy, openpyxl, matplotlib가 설치된 Python 환경을 사용한다.

```bash
python -m sources.mcm_perturbation.score_samples
```

저장된 점수로 그림 24–25를 다시 그릴 때는 아래 명령을 사용한다.

```bash
python -m sources.mcm_perturbation.plot_program_scores
```

로컬 발현은 `data/perturbation_validation/`에서 읽고 GEO metadata 및 orthology는 `data/perturbation_validation/reference/`의 작은 참조 파일을 사용한다. 이 참조의 원본과 해시는 [provenance.tsv](../../data/perturbation_validation/reference/provenance.tsv)에 있다. Mouse에는 로컬 `datamining_5`에서 검증한 Ensembl Compara 후보표의 high-confidence one-to-one 관계를 사용하고, stable ID는 human GRCh38.114 및 mouse GRCm39.116 GTF로 symbol에 연결했다. 두 표준 annotation에 각각 유일한 symbol만 연결했다. ACHES의 `HLA.E`는 `HLA-E`로 연결했다. SNAP loading의 출처는 Ling et al. Supplementary Tables 6–7이며 공개 loading의 부호를 유지했다. [Ling et al.](https://doi.org/10.1038/s41586-024-07109-5)

각 코호트에서 예정 비교 표본 전체에서 0인 유전자를 제외했다. GSE228927에는 log2(FPKM+1), GSE282665에는 log2(TPM+1)을 적용했다. ACHES는 유전자 평균, SNAP-a와 SNAP-n은 signed weighted sum을 retained absolute loading 합으로 나눈 값이다.

| 결과 | 내용 |
| --- | --- |
| [planned_sample_scores.tsv](../../data/perturbation_validation/results/planned_sample_scores.tsv) | 예정 비교에 포함된 표본별 점수와 GSM accession |
| [program_scores.tsv](../../data/perturbation_validation/results/program_scores.tsv) | Human recipient matrix의 모든 열과 예정 mouse MCM 열의 점수 및 원래 가중합 |
| [sample_manifest.tsv](../../data/perturbation_validation/results/sample_manifest.tsv) | 원본 발현 열의 GEO 대응, 조건, recipient와 source 정보, 포함 여부 |
| [paired_score_differences.tsv](../../data/perturbation_validation/results/paired_score_differences.tsv) | Human recipient line별 LPS-MCM − resting-MCM 차이 |
| [condition_summary.tsv](../../data/perturbation_validation/results/condition_summary.tsv) | 예정 비교의 조건별 표본 열 수, 평균, SD |
| [program_coverage.tsv](../../data/perturbation_validation/results/program_coverage.tsv) / [program_mapping.tsv.gz](../../data/perturbation_validation/results/program_mapping.tsv.gz) | 사용 유전자 수, 전체 및 부호별 가중치 coverage, 유전자별 매핑 상태 |
| [program_score_effects_ci.tsv](../../data/perturbation_validation/results/program_score_effects_ci.tsv) | 처리군 − resting-MCM의 평균 차이, 표준오차, 자유도, 95% 신뢰구간, p값과 q값 |
| GSE282665 괄호 그림 · [PDF](../../figures/24.pdf) · [SVG](../../figures/24.svg) · [PNG](../../figures/24.png) | Human astrocyte의 ACHES와 SNAP-a 표본 점수 및 line 내 LPS 대비 |
| GSE228927 괄호 그림 · [PDF](../../figures/25.pdf) · [SVG](../../figures/25.svg) · [PNG](../../figures/25.png) | Mouse astrocyte의 ACHES와 SNAP-a 표본 점수 및 resting-MCM 대비 LPS-MCM과 hemin-MCM |

Human 신뢰구간은 recipient line 내 차이에 paired t 방식을 적용해 계산한다. Motor neuron의 기계산 점수와 대비는 전체 TSV에 보관한다. [GSE282665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282665)

GSE228927의 `NOT` 표본은 resting-MCM 대조군이다. LPS와 hemin은 각각 같은 resting-MCM 표본과 비교했으며, 조건 간 pairing은 확인되지 않아 독립군을 가정한 양측 Welch t 검정을 사용했다. Mouse SNAP-a는 원본 유전자의 41%와 절대 가중치의 63%를 사용했으므로 탐색적 결과로 표시했다. [GSE228927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228927)

두 코호트의 astrocyte LPS 대비 네 개를 함께 BH 보정했다. Mouse hemin 대비 두 개와 human motor neuron SNAP-n 대비는 별도로 보정했다. `*`는 q < 0.05, `**`는 q < 0.01, `***`는 q < 0.001, `ns`는 q ≥ 0.05다.
