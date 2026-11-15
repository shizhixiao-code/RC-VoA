"""Standalone experiment report; never writes manuscript or Figure 4 files."""
from pathlib import Path
import hashlib
import html
import json
import numpy as np
import pandas as pd
import scipy.optimize._highspy._core as highs
from run_experiment import load, exact_args
from rc_interval import certificate
from exact_completion import NAMES

HERE = Path(__file__).resolve().parent


def main():
    df = pd.read_csv(HERE/'full_results.csv')
    online = pd.read_csv(HERE/'online_results.csv')
    assert len(df) == 20000 and df.table_id.nunique() == 200
    assert not df[['table_id','order']].duplicated().any()
    assert df.groupby('table_id').size().eq(100).all()
    assert df.verdict.eq(df.truth).all() and df.solver_unknowns.eq(0).all()
    assert df.exact_stop.le(df.rc_stop).all() and df.earliest_proved.all()
    checked = df[~df.screen].copy()
    assert checked.table_id.nunique() == 56
    online = online.merge(df[['table_id','order','exact_stop','rc_stop','exact_checks','rc_checks']],
                          on=['table_id','order'],validate='many_to_one')
    for method, stop, checks in [('Exact MILP','exact_stop','exact_checks'),('RC-VoA','rc_stop','rc_checks')]:
        part = online[online.method.eq(method)]
        proved = part[part.unknowns.eq(0)]
        assert proved.stop.eq(proved[stop]).all() and proved.checks.eq(proved[checks]).all()
        assert part.stop.ge(part[stop]).all() and part.checks.ge(part[checks]).all()
    order_rows, summaries, detail_rows = [], [], []
    for stage, feasible in [('Primary',False),('Primary',True),('Repaired',False),('Repaired',True)]:
        part = checked[checked.stage.eq(stage)&checked.feasible.eq(feasible)]
        po = part.groupby('order')[['rc_fraction','exact_fraction']].mean().mul(100)
        po['gain_pp'] = po.rc_fraction-po.exact_fraction
        for oi, r in po.iterrows():
            order_rows.append(dict(stage=stage, feasible=feasible, order=oi, **r.to_dict()))
        timing = online[online.stage.eq(stage)&online.feasible.eq(feasible)]
        row = dict(stage=stage, feasible=feasible, configurations=part.table_id.nunique(),
            comparisons=len(part), rc_mean_pct=po.rc_fraction.mean(), rc_p10=po.rc_fraction.quantile(.1),
            rc_p90=po.rc_fraction.quantile(.9), exact_mean_pct=po.exact_fraction.mean(),
            exact_p10=po.exact_fraction.quantile(.1), exact_p90=po.exact_fraction.quantile(.9),
            full_mean_pct=100.0, gain_pp=po.gain_pp.mean(), gain_p10=po.gain_pp.quantile(.1),
            gain_p90=po.gain_pp.quantile(.9), mean_saved_records=part.record_saving.mean(),
            earlier_orders=int(part.exact_stop.lt(part.rc_stop).sum()),
            same_orders=int(part.exact_stop.eq(part.rc_stop).sum()), max_saved_records=int(part.record_saving.max()),
            max_gain_pp=float((100*(part.rc_fraction-part.exact_fraction)).max()),
            max_saved_groups=int((part.rc_stop-part.exact_stop).max()))
        for method, short in [('RC-VoA','rc'),('Exact MILP','exact'),('Pre-screen + full','full')]:
            t = timing[timing.method.eq(method)]
            row[short+'_online_n'] = len(t)
            row[short+'_online_mean_s'] = t.wall_seconds.mean()
            row[short+'_online_p10_s'] = t.wall_seconds.quantile(.1)
            row[short+'_online_p90_s'] = t.wall_seconds.quantile(.9)
            row[short+'_online_max_s'] = t.wall_seconds.max()
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary.to_csv(HERE/'summary.csv',index=False)
    pd.DataFrame(order_rows).to_csv(HERE/'summary_by_order.csv',index=False)
    config = checked.groupby(['table_id','stage','cohort','mode','batch','cap','feasible'],as_index=False).agg(
        rc_fraction=('rc_fraction','mean'),exact_fraction=('exact_fraction','mean'),
        mean_saved_records=('record_saving','mean'),max_saved_records=('record_saving','max'),
        earlier_orders=('record_saving',lambda x:int(x.gt(0).sum())))
    config['gain_pp'] = 100*(config.rc_fraction-config.exact_fraction)
    config.to_csv(HERE/'configuration_summary.csv',index=False)
    negative_stops = []
    for key, g in checked[~checked.feasible].groupby('table_id'):
        data = load(key)
        for rec in g.itertuples():
            seen = data['group_counts'][data['orders'][rec.order, :rec.exact_stop]].sum(0)
            v, lo, hi = certificate(*exact_args(data, seen))
            negative_stops.append(dict(table_id=key, stage=rec.stage, order=rec.order,
                                       **{name:bool(hi[j]<0) for j,name in enumerate(NAMES)}))
    negatives = pd.DataFrame(negative_stops)
    negatives.to_csv(HERE/'infeasible_stop_reasons.csv',index=False)
    # Largest absolute endpoint saving is a labeled descriptive example, not a
    # cherry-picked replacement for the complete configuration/order summaries.
    best = checked.loc[checked.record_saving.idxmax()]
    data = load(best.table_id)
    seen = data['group_counts'][data['orders'][int(best.order), :int(best.exact_stop)]].sum(0)
    v, lo, hi = certificate(*exact_args(data,seen))
    example = dict(table_id=best.table_id, order=int(best.order), eligible=int(best.eligible),
        exact_prefix=int(best.exact_stop), rc_prefix=int(best.rc_stop), exact_checks=int(best.exact_checks),
        rc_checks=int(best.rc_checks), saved_records=int(best.record_saving),
        remaining_eligible=int(best.eligible-best.exact_checks), rc_verdict_at_exact_stop=int(v),
        rc_lower_margins=dict(zip(NAMES,lo.tolist())), rc_upper_margins=dict(zip(NAMES,hi.tolist())))
    (HERE/'largest_saving_example.json').write_text(json.dumps(example,indent=2))
    # Compare exact feasible costs with the archived necessary-record lower bound.
    source = HERE.parents[1]/'rcvoa_threshold_transport_20261102/analysis/rescue_verification_extension_20261026/verification'
    lb = pd.read_csv(source/'positive_1400_bounds.csv')
    bound = checked[checked.stage.eq('Repaired')&checked.feasible].merge(lb[['table_id','order','lower_bound']],
                    on=['table_id','order'],validate='one_to_one')
    assert len(bound)==1400 and bound.exact_checks.ge(bound.lower_bound).all()
    bounds = dict(repaired_necessary_mean_pct=float((100*bound.lower_bound/bound.eligible).mean()),
                  repaired_exact_mean_pct=float((100*bound.exact_fraction).mean()),
                  necessary_bound_attained_orders=int(bound.exact_checks.eq(bound.lower_bound).sum()),
                  necessary_bound_total_orders=len(bound))
    validation = json.loads((HERE/'VALIDATION.json').read_text())
    audit = json.loads((HERE/'BOUNDARY_AUDIT.json').read_text())
    assert validation['status']==audit['status']=='PASS'
    solver = '.'.join(str(getattr(highs, 'HIGHS_VERSION_'+x)) for x in ['MAJOR','MINOR','PATCH'])
    bench = json.loads((HERE/'BENCHMARK_ENVIRONMENT.json').read_text(encoding='utf-8-sig'))
    summary_json = dict(rows=len(df), nonscreen_rows=len(checked), agreements=int(df.verdict.eq(df.truth).sum()),
        solver_unknowns=int(df.solver_unknowns.sum()+online.unknowns.sum()),
        first_stop_equal=int(df.exact_stop.eq(df.rc_stop).sum()), first_stop_earlier=int(df.exact_stop.lt(df.rc_stop).sum()),
        highs_version=solver, benchmark_environment=bench, online_pairs=len(online)//3, online_orders=sorted(online.order.unique().tolist()),
        strata=summary.to_dict('records'), bounds=bounds, validation=validation, independent_audit=audit)
    summary_json['infeasible_stopping_reasons_nonexclusive'] = negatives.groupby('stage')[NAMES].sum().to_dict('index')
    summary_json['largest_absolute_saving_example'] = example
    (HERE/'RESULTS.json').write_text(json.dumps(summary_json,indent=2),encoding='utf-8')
    title = 'Exact-Completion MILP：独立实验结果'
    lines = [title, '日期：2026-11-07（Asia/Shanghai）',
        '范围：仅实验；论文、原 Figure 4、原始输入均未修改。', '',
        '主要发现：需检索的4,100次infeasible判定，两种方法的最早停止组完全一致。',
        'Repaired feasible平均查询比例99.984012%→99.626972%，节省0.357040个百分点（5.472条记录/次）；',
        '419/1,400个顺序提前停止，最大单次节省103条（6.706个百分点、8个组）。',
        '同层在线均值0.023571秒→1.629214秒；计算更慢，但绝对耗时是秒级，实际取舍取决于诊断检索成本。',
        '这些结果仅说明所测固定顺序下的保守性，不证明全局最优查询策略。', '',
        '完整配对：100 primary + 100 repaired，各100个原有 outcome-blind group orders，共20,000次。',
        '20% imaging cap、free quota、同一 complete criterion、已知 eligibility、无 outcome totals。',
        '两阶段各72个同一配置经共同 decision-only pre-screen 以零查询排除，以下不混入这72个配置。', '',
        '查询比例：先在各表内除以 eligible endpoints，再等权平均；括号为跨100个配对顺序的P10–P90（非置信区间）。']
    for r in summaries:
        tag = r['stage']+' '+('feasible' if r['feasible'] else 'infeasible')
        lines.append(f"{tag} (n={r['configurations']}): RC {r['rc_mean_pct']:.6f}% ({r['rc_p10']:.6f}–{r['rc_p90']:.6f}); "
                     f"Exact {r['exact_mean_pct']:.6f}% ({r['exact_p10']:.6f}–{r['exact_p90']:.6f}); Full 100%; "
                     f"gain {r['gain_pp']:.6f} percentage points; {r['earlier_orders']}/{r['comparisons']} orders earlier; "
                     f"mean {r['mean_saved_records']:.4f} fewer endpoints.")
    lines += ['', '实际逐组在线计算时间（不含真实病理检索延迟）：']
    for r in summaries:
        tag = r['stage']+' '+('feasible' if r['feasible'] else 'infeasible')
        lines.append(f"{tag}: {r['exact_online_n']} table–order runs; mean RC {r['rc_online_mean_s']:.6f}s, "
                     f"Exact {r['exact_online_mean_s']:.6f}s, pre-screen+full {r['full_online_mean_s']:.6f}s; "
                     f"Exact run P10–P90 {r['exact_online_p10_s']:.6f}–{r['exact_online_p90_s']:.6f}s.")
    lines += [f"在线采样预设为全部56个需查询配置的原始顺序 {summary_json['online_orders']}；未用二分搜索或跨前缀缓存。",
              '在线时间包含 ledger 累加、模型构造和求解；不含磁盘载入、Python 启动和真实诊断获取。',
              'Pre-screen + full 的时间只测本地模拟累加和最后一次计数判定，因此不代表实际检索服务时间。',
              '不同方法在同一进程逐一运行；当前通用实现未经独立性能调优，时间依赖本机和求解器。', '',
              f"一致性：{summary_json['agreements']}/20,000；本次运行 UNKNOWN/timeout = {summary_json['solver_unknowns']}。",
              f"独立验证：{validation['completions']}种补全、{validation['individual_margin_checks']}个逐margin检查、"
              f"{validation['literal_allocation_completions']}个literal allocation穷举。",
              f"最早性证据：{audit['opposite_completion_pairs']}个 stop-1 前缀均保存一份可行补全和一份不可行补全；"
              f"独立分配MILP共{audit['independent_milps']}次，全部一致。",
              f"Repaired feasible必要下界：{bounds['repaired_necessary_mean_pct']:.6f}%；Exact：{bounds['repaired_exact_mean_pct']:.6f}%；"
              f"{bounds['necessary_bound_attained_orders']}/1400达到该下界。", '',
              '实现与适用边界：',
              '1. x[j,c]为7类合法诊断、8个(B,Dpair,R)组合的非负整数计数，x>=seen且每列总和固定为seen+unresolved。',
              '   同一组合内记录的系数完全相同，任一整数计数都能分配给未查询的个体，故聚合不缩小补全空间。',
              '   group只限制哪些记录一起揭示，不要求同组记录诊断相同。Excluded不是可选的未知诊断。',
              '2. 可行补全采用全部margin同时成立的整数epigraph模型；将deficit/max变量取最小精确值只会改善各margin，故存在性等价。',
              '3. 不可行补全逐margin搜索；ceil由0 <= d*y-num <= d-1精确表示，max使用有界binary disjunction，整数负margin写为<=-1。',
              '4. 求解器仅status=2视为UNSAT，status=0的可行解还需逐项整数重验；time/node limit、错误或未知均为UNRESOLVED。',
              '   本次单个MILP时限10s，gap=0；binary offline search遇未知会中止而不声称最早性。在线模式遇未知继续检索。',
              '5. 停止点采用离线二分查找并验证 stop-1 的双补全证据；离线缓存/二分耗时不作为部署时间。',
              '6. 只覆盖Figure 4的at-most-cap/free-quota，不自动推广到固定非零quota或exact-volume扩展；也不声称全局最优查询顺序。',
              '7. 配置共享病例、模型和校准，不是独立临床样本；未对order variation做显著性/临床收益推断。', '',
              '运行环境：Python 3.13.14，NumPy 2.5.1，SciPy 1.18.0，HiGHS '+solver+'；CPU：'+bench['cpu']+'。',
              '一个Python进程，方法依次执行；HiGHS线程数沿用默认值，未强制限制为单线程。',
              'SciPy状态定义：https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.milp.html',
              '复现顺序：prepare.py → validate.py → run_experiment.py pilot → run_experiment.py full →',
              '          audit_boundaries.py → run_experiment.py online --orders 10 → summarize.py',
              'prepare.py需要工作区原归档；其余脚本使用本目录已封存inputs。',
              '详情：summary.csv、configuration_summary.csv、full_results.csv、online_results.csv、offline/*_boundary_witnesses.json。']
    (HERE/'RESULTS_ZH.txt').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    table_rows = []
    for r in summaries:
        tag=r['stage']+' '+('feasible' if r['feasible'] else 'infeasible')
        table_rows.append('<tr>'+''.join(f'<td>{c}</td>' for c in [tag,str(r['configurations']),
            f"{r['rc_mean_pct']:.4f}%<small>{r['rc_p10']:.4f}–{r['rc_p90']:.4f}</small>",
            f"{r['exact_mean_pct']:.4f}%<small>{r['exact_p10']:.4f}–{r['exact_p90']:.4f}</small>",
            '100%',f"{r['gain_pp']:.4f} pp", f"{r['earlier_orders']}/{r['comparisons']}"] )+'</tr>')
    timing_rows=[]
    for r in summaries:
        timing_rows.append('<tr>'+''.join(f'<td>{c}</td>' for c in [r['stage']+' '+('feasible' if r['feasible'] else 'infeasible'),
            str(r['exact_online_n']),f"{1000*r['rc_online_mean_s']:.2f} ms",f"{1000*r['exact_online_mean_s']:.2f} ms",
            f"{1000*r['exact_online_p10_s']:.2f}–{1000*r['exact_online_p90_s']:.2f} ms"] )+'</tr>')
    document = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>'+title+'</title><style>'
    document += 'body{font:16px/1.65 system-ui,"Microsoft YaHei",sans-serif;color:#172c3b;background:#f6f8fa;max-width:1180px;margin:40px auto;padding:0 24px}h1{font-size:30px}h2{font-size:22px;margin-top:32px}table{border-collapse:collapse;width:100%;background:white}th,td{text-align:left;padding:12px;border-bottom:1px solid #dce3e9}th{background:#eaf0f5}small{display:block;color:#586a78}pre{white-space:pre-wrap;background:white;padding:24px;border:1px solid #dce3e9}a{color:#176ab6}.badge{padding:12px 16px;border-left:4px solid #198164;background:#e9f5ef}</style>'
    document += f'<h1>{title}</h1><p>2026-11-07 · 同一 Figure 4 条件 · 论文未修改</p><div class="badge">20,000/20,000 判定一致 · {summary_json["solver_unknowns"]} 次超时/未知 · 56 个需查询配置完整配对 · 72 个零查询配置单独统计</div>'
    document += '<p>需检索的4,100次不可行判定全部同组停止。改进主要出现在可行配置：repaired feasible平均少查5.47条记录（0.3570个百分点），419/1,400个顺序提前停止，最大单次少查103条。其在线计算均值由0.0236秒增至1.629秒；计算与检索的实际取舍仍取决于记录获取成本。</p>'
    document += f'<p>独立复核：{audit["opposite_completion_pairs"]:,}个停止前缀的双补全证据，{audit["independent_milps"]:,}次独立MILP检查，全部通过。Repaired feasible精确均值99.6270%，仍高于99.3051%的必要下界；1,254/1,400个顺序达到该下界。</p>'
    document += '<h2>需要检索的配置</h2><p>均值和括号内 P10–P90 均基于先归一化、再等权平均的100个配对查询顺序；P10–P90不是置信区间。</p>'
    document += '<table><thead><tr>'+''.join('<th>'+h+'</th>' for h in ['分层','n','RC-VoA','Exact MILP','Pre-screen + full','额外节省','更早停止'])+'</tr></thead><tbody>'+''.join(table_rows)+'</tbody></table>'
    document += '<h2>真正逐组运行的计算时间</h2><p>全部56个需查询配置，各取预设前10个顺序；包括模型构造与求解，不含实际病理检索时间。</p>'
    document += '<table><thead><tr>'+''.join('<th>'+h+'</th>' for h in ['分层','运行数','RC均值','Exact均值','Exact运行P10–P90'])+'</tr></thead><tbody>'+''.join(timing_rows)+'</tbody></table>'
    document += '<p><a href="summary.csv">分层CSV</a> · <a href="full_results.csv">20,000条配对明细</a> · <a href="configuration_summary.csv">各配置结果</a> · <a href="online_results.csv">在线计时明细</a> · <a href="RESULTS.json">机器可读报告</a></p>'
    document += '<h2>方法、验证和解释边界</h2><pre>'+html.escape('\n'.join(lines[lines.index('实现与适用边界：'):]))+'</pre></html>'
    (HERE/'RESULTS_ZH.html').write_text(document,encoding='utf-8')
    manifest = {str(p.relative_to(HERE)):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.json'}
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(summary.to_string(index=False),flush=True)
    print(json.dumps(bounds,indent=2),flush=True)


if __name__=='__main__': main()
