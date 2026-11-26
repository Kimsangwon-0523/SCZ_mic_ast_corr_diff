# MCM perturbation 분석

현재 그림과 결과 노트북은 GSE228927과 GSE282665의 astrocyte ACHES 및 SNAP-a를 보여준다. [그림 색인](../../FIGURES.md)에서 번호 24–25을 확인할 수 있다.

- [00_inspect_perturbation_data.ipynb](00_inspect_perturbation_data.ipynb)에서 발현 파일과 원본 열을 확인한다.
- [01_review_program_scores.ipynb](01_review_program_scores.ipynb)에서 데이터셋별 표본 점수, 비교 괄호의 별표와 q값, 조건 간 차이의 95% 신뢰구간을 확인한다.

## 다시 계산하기

저장소 루트에서 Python에 pandas, numpy, scipy, openpyxl, matplotlib가 있고 R에 edgeR가 설치된 환경을 사용한다.

```bash
python -m sources.mcm_perturbation.score_samples --rscript /path/to/Rscript
```

저장된 점수로 두 데이터셋의 그림만 다시 그릴 때는 아래 명령을 사용한다.

```bash
python -m sources.mcm_perturbation.plot_program_scores
```

로컬 발현은 `data/perturbation_validation/`에서 읽고 GEO metadata 및 orthology는 `data/perturbation_validation/reference/`의 작은 참조 파일을 사용한다. 이 참조의 원본과 해시는 [provenance.tsv](../../data/perturbation_validation/reference/provenance.tsv)에 있다. Mouse에는 로컬 `datamining_5`에서 검증한 Ensembl Compara 후보표의 high-confidence one-to-one 관계를 사용하고, stable ID는 human GRCh38.114 및 mouse GRCm39.116 GTF로 symbol에 연결했다. 두 표준 annotation에 각각 유일한 symbol만 연결했다. ACHES의 `HLA.E`는 `HLA-E`로 연결했다. SNAP loading의 출처는 Ling et al. Supplementary Tables 6–7이며 공개 loading의 부호를 유지했다. [Ling et al.](https://doi.org/10.1038/s41586-024-07109-5)

GSE175485에서는 MCM 6열에 edgeR `filterByExpr`, TMM, logCPM(prior count 1)을 적용했다. GSE228927은 log2(FPKM+1), GSE282665는 log2(TPM+1)이다. 후자의 두 데이터는 recipient별로 예정 비교 표본 전체에서 0인 유전자를 제외한다. 그때 확정한 gene set을 같은 receiver matrix의 모든 점수에 적용했다. ACHES는 평균, SNAP-a와 SNAP-n은 signed weighted sum을 retained absolute loading 합으로 나눈 값이다. Coverage가 기준보다 낮으면 동일한 계산 결과를 기술 요약과 함께 표시한다.

| 결과 | 내용 |
| --- | --- |
| [planned_sample_scores.tsv](../../data/perturbation_validation/results/planned_sample_scores.tsv) | 예정 비교에 포함된 표본별 점수와 GSM accession |
| [program_scores.tsv](../../data/perturbation_validation/results/program_scores.tsv) | Human recipient matrix의 모든 열과 예정 mouse MCM 열의 점수 및 원래 가중합 |
| [sample_manifest.tsv](../../data/perturbation_validation/results/sample_manifest.tsv) | 원본 발현 열의 GEO 대응, 조건, recipient와 source 정보, 포함 여부 |
| [paired_score_differences.tsv](../../data/perturbation_validation/results/paired_score_differences.tsv) | Human recipient line별 LPS-MCM − resting-MCM 차이 |
| [condition_summary.tsv](../../data/perturbation_validation/results/condition_summary.tsv) | 예정 비교의 조건별 표본 열 수, 평균, SD |
| [program_coverage.tsv](../../data/perturbation_validation/results/program_coverage.tsv) / [program_mapping.tsv.gz](../../data/perturbation_validation/results/program_mapping.tsv.gz) | 사용 유전자 수, 전체 및 부호별 가중치 coverage, 유전자별 매핑 상태 |
| [program_score_effects_ci.tsv](../../data/perturbation_validation/results/program_score_effects_ci.tsv) | 처리군 − resting-MCM의 평균 차이, 표준오차, 자유도, 95% 신뢰구간, p값과 q값 |
| GSE228927 괄호 그림 · [PDF](../../figures/24.pdf) · [SVG](../../figures/24.svg) · [PNG](../../figures/24.png) | Mouse astrocyte ACHES와 SNAP-a의 표본 점수, LPS와 hemin 대비 |
| GSE282665 괄호 그림 · [PDF](../../figures/25.pdf) · [SVG](../../figures/25.svg) · [PNG](../../figures/25.png) | Human astrocyte의 ACHES와 SNAP-a 표본 점수 및 line 내 LPS 대비 |

Human 신뢰구간은 recipient line 내 차이에 paired t 방식을 적용해 계산한다. Motor neuron의 기계산 점수와 대비는 전체 TSV에 보관한다. [GSE282665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282665)

Mouse GEO 기록에는 별도 실험 3회가 명시되지만 조건 간 실험 매칭은 제공되지 않는다. 각 대비의 Welch 신뢰구간은 실험군과 대조군의 독립을 가정한 탐색적 계산이며 두 대비가 resting-MCM 표본을 공유한다. 신뢰구간은 대비별 95% 구간이다. [GSE228927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228927)

괄호의 별표는 양측 t 검정의 BH 보정 q값을 사용한다. 현재 두 데이터셋의 astrocyte LPS 비교 4개를 primary family로, mouse hemin 2개와 human motor neuron LPS 1개를 secondary family로 보정했다. `*`는 q < 0.05, `**`는 q < 0.01, `***`는 q < 0.001, `ns`는 q ≥ 0.05다. Mouse SNAP-a의 낮은 유전자 coverage는 그림에 †로 표시했다.
