# 阈值迁移实验：完整结果与证据边界

原始 source-S/P、36 张主表、1→14 修复、所有既有证明及核查结果保留。新增实验不训练模型，不替换 source-S。

## 原始阈值与数据核验

- 24 个 source S/P 阈值均由匹配 MILK 分数逐值复现；连续分数覆盖全部 9,168 个外部病灶。
- 9 个配对系统使用 3 个不同临床分数向量；L1–L3 不能算作独立临床阈值证据。
- S 的 MILK 整数要求：melanoma 419/441、SCC 296/311、BCC 2259/2509、IEC 162/179；另需 pooled 3105/3449。P 只要求 pooled 3105/3449。三条临床 S 阈值均由 melanoma 约束决定。
- 原有 36 份分组划分彼此一致；所有 calibration/evaluation groups 均不相交。MRA 使用原记录中的 groups，其他队列使用已归档 patient/participant 分组。

## Oracle 与独立本地校准

- 满足全部观察到的敏感度目标时，UQ 单阈值 oracle 标记率仍为 98.54%–99.45%；MRA 为 92.70%–96.86%。
- Oracle 使用 evaluation labels，仅界定固定分数、单阈值类的回顾性最佳工作点。它没有参与 S_local 选择。
- S_local 的 60 个阈值先由 calibration folds 选择并冻结，之后才评估。UQ 的 melanoma 有 16 个病灶、14 个参与者；95% 整数目标要求检出全部 16 个。MCR-SL 的 8 个 melanoma 也必须全部检出。
- ConvNeXt/UQ：标记率从 99.64% 降至 70.67%，FP 从 7552 降至 5336；pooled 由 95/95 降至 93/95，melanoma 由 16/16 降至 14/16。该本地阈值没有满足全部敏感度目标。
- ConvNeXt/MRA：标记率从 99.78% 降至 92.93%，FP 从 422 降至 373；pooled 377/379、melanoma 68/70，聚合敏感度目标达标，但丢失的两个 source-S 癌症仍必须被保护。
- 临床分数不是完全没有区分力：UQ 的 pooled AUROC 分别为 0.673、0.837、0.814；但总体排序能力不保证严格亚型敏感度尾部仍有足够特异性。本实验区分了阈值调整空间与固定分数的受约束工作量，不能独立识别分布偏移的因果来源。

## 不确定性与稀有亚型

- 1,000 次原分组 bootstrap，按 evaluation fold 分层；source/local 使用相同抽样。区间以模型和阈值已固定为条件，不包括阈值估计不确定性。
- ConvNeXt/UQ 的标记率差为 −28.97 pp [−31.91, −26.02]；local melanoma 敏感度 87.5% [66.7, 100.0]。达标的 calibration fold 不提供 evaluation 或人群敏感度保证。
- 每个亚型的病例数、分组数、无该亚型的 bootstrap 次数及退化区间均单列。不存在于整个队列的亚型与存在于队列但缺乏 calibration 证据的亚型明确区分。

## 保留 source-S 保护后的完整准则

- 令 B=S_local、D_pair=原始 source-P、R=source-S。B/R 不需要相等或嵌套；固定阈值保持 selection invariance。六类裕度沿用双基线扩展，不能只套用单基线固定目标命题。
- 36 张合并 cross-fitted 表在 15%/20% 均为 0/36 可行；每张均缺乏至少一例 mandatory rescue。原始 source-S/P 同预算仍为 27/36 可行。两种 ranker 在合并表上都没有通过。
- 同记录逐折对照：source-S 为 165/180（15%/20%）；S_local 为 96/180（15%）、99/180（20%）。20% 下 local Proximity/Clinical-risk 通过 12/180、8/180，source 分别为 46/180、17/180。两种审计粒度不能合并或互相替换。
- Local 的非零 rescue 正例来自 UQ 第 5 个 evaluation fold 的 DeiT L1–L3：t*=1、C_minus=1、m=105。这是 3 个共享临床分数的配置、1 个患者组，不是 3 个独立临床复制，也不能加到原来的 14/100 分母中。
- 3,327 次独立逐病灶整数规划与公式一致；9,072 个区间前缀、1,944 个替代完成表及 1,296 个完整表/预算检查均通过。既有核查成本结论和 23,600 的原分母不变。

## 完整记录

- `Table_SIV_full.csv`：三种设置全部 36 个独特临床行，所有 TP/n、FP、flag、阈值及每亚型结果。
- `all_36_system_operating_points.csv`：到 9 个系统 × 4 队列的完整映射，明确共享临床分数。
- `operating_points.csv` / `calibration_coverage.csv` / `uncertainty.csv`：逐折评估、校准缺口及成对区间。
- `complete_criterion_tables.csv` / `subtype_rescue_margins.csv` / `ranker_results.csv`：完整准则、六类裕度、合法 quota 与实际策略。
- `minimum_witnesses.csv.gz` / `policy_selections.csv.gz`：病灶级 witness 与选择记录。
- `frontier.csv`：完整 whole-tie 枚举，包含未出现在印刷表中的所有工作点。
- Preparation seal、源文件 hash、数学适用性说明、MILP 和区间检查位于实验目录。
