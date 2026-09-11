#!/usr/bin/env python3
"""Export paper Tables 2--4 from the reference or regenerated records.

All means and binomial intervals are recomputed from per-stream data. This
presentation layer does not alter scientific results or sampling protocols.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from scipy.stats import beta


def cp(count: int, total: int) -> tuple[float,float]:
    return (float(beta.ppf(.025,count,total-count+1)) if count else 0.,
            float(beta.ppf(.975,count+1,total-count)) if count<total else 1.)


def markdown(headers: list[str], rows: list[list[str]]) -> str:
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(row)+' |' for row in rows])


def latex(headers: list[str], rows: list[list[str]], align: str) -> str:
    return '\n'.join([r'\begin{tabular}{'+align+'}', r'\toprule', ' & '.join(headers)+r'\\',r'\midrule']+
                     [' & '.join(row)+r'\\' for row in rows]+[r'\bottomrule',r'\end{tabular}',''])


def export(workspace: Path) -> None:
    out=workspace/'tables';out.mkdir(exist_ok=True)
    endpoints=pd.read_csv(workspace/'data/horizon_endpoint_comparison.csv')
    table2=endpoints.loc[endpoints.N==100000,['k','beta','CG23','finite_horizon','anytime','any_union']].copy()
    table2.columns=['k','n_Beta','n_CG23','n_FH','n_Anytime','Anytime_over_FH']
    table2.to_csv(out/'table2.csv',index=False)
    rows2=[[str(int(r.k)),f'{r.n_Beta:.2f}',f'{r.n_CG23:.2f}',f'{r.n_FH:.2f}',f'{r.n_Anytime:.2f}',f'{r.Anytime_over_FH:.3f}'] for r in table2.itertuples()]
    headers2=['$k$','Beta','CG23','FH','Any','Any/FH']
    (out/'table2.tex').write_text(latex(headers2,rows2,'rrrrrr'))
    first=pd.read_csv(workspace/'data/stopping_audit_replicates.csv')
    horizon=pd.read_csv(workspace/'data/horizon_stopping_replicates.csv')
    specs=[('Fixed N','Fixed N',first),('CG23-repeated','CG23-repeated',first),
           ('GCC--T4','GCC refined',first),('FH, N=800','Finite horizon',horizon),
           ('Summable CG23','Summable CG23',horizon),('Anytime','Anytime',horizon)]
    records=[];rows3=[]
    for label,name,df in specs:
        selected=df[(df.method==name)&(df.eps==.05)]
        if len(selected)!=10000 or int(selected.deployed.sum())!=10000:
            raise ValueError(f'{label}: expected all 10,000 streams to return a design.')
        errors=int(selected.risk_exceeds_eps.sum());lo,hi=cp(errors,len(selected))
        mean=float(selected.acquired.mean())
        records.append(dict(method=label,streams=len(selected),mean_n=mean,errors=errors,error_probability=errors/len(selected),CP95_lower=lo,CP95_upper=hi))
        if errors==0:interval=f'0 [0, {100*hi:.4f}]'
        elif errors>=100:interval=f'{errors/100:.2f} [{100*lo:.2f}, {100*hi:.2f}]'
        else:interval=f'{errors/100:.2f} [{100*lo:.3f}, {100*hi:.3f}]'
        rows3.append([label,f'{mean:.1f}',str(errors),interval])
    pd.DataFrame(records).to_csv(out/'table3.csv',index=False)
    headers3=['Method','Mean $n$','Errors',r'Error probability (\%)']
    (out/'table3.tex').write_text(latex(headers3,rows3,'lrrr'))
    monitoring=pd.read_csv(workspace/'data/horizon_monitoring_8000.csv')
    if len(monitoring)!=10000:raise ValueError('Expected the full 10,000-stream monitoring cohort.')
    records4=[];rows4=[]
    for name,label in [('CG23','CG23'),('Finite horizon','FH, N=8000'),('Summable CG23','Summable CG23'),('Anytime','Anytime')]:
        row=[label]
        for n in (800,8000):
            first_crossing=monitoring[name]
            count=int(((first_crossing>0)&(first_crossing<=n)).sum())
            lo,hi=cp(count,len(monitoring))
            records4.append(dict(method=label,horizon=n,streams=len(monitoring),crossings=count,crossing_probability=count/len(monitoring),CP95_lower=lo,CP95_upper=hi))
            row.append(str(count))
        rows4.append(row)
    pd.DataFrame(records4).to_csv(out/'table4.csv',index=False)
    headers4=['Bound','Through 800','Through 8000']
    (out/'table4.tex').write_text(latex(headers4,rows4,'lrr'))
    sections=[('# Numerical results\n\nThis package generates numerical results for the paper '
               '**A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design**, by **G. Calafiore**.\n\n'
               'Tables below are recomputed from this workspace. '
               'The numerical verification scripts separately check the printed values.\n'),
              '## Table 2: certificate comparisons\n\n'+markdown(headers2,rows2),
              '## Table 3: first-target stopping\n\n'+markdown(['Method','Mean n','Errors','Error probability (%)'],rows3),
              '## Table 4: continued inspection\n\n'+markdown(headers4,rows4)]
    for name,label in [('horizon_comparison','Figure 1: horizon comparison'),('certificate_archive','Figure 2: certified archive'),('disturbance_risk','Figure 3: disturbance geometry')]:
        if (workspace/'figures'/f'{name}.png').is_file():
            sections.append(f'## {label}\n\n![{label}](figures/{name}.png)')
    sections.append('## Interpretation\n\nThe mechanical candidate bound remains conditional on almost-sure support reconstruction. '
                    'Its validation counts contraction/LMI violations. The aligned and alternative performance evaluators retain their source definitions; '
                    'the alternative comparison is exploratory. These runs use the original certificate, without the persistence-window minimum.')
    (workspace/'paper_results.md').write_text('\n\n'.join(sections)+'\n')
    print('Exported tables/table{2,3,4}.{csv,tex} and paper_results.md.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',type=Path,required=True)
    args=parser.parse_args();export(args.workspace.resolve())

if __name__=='__main__':main()
