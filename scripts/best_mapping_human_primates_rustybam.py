import numpy as np
import pandas as pd
import argparse
import re

def get_length(pos):
    pos = pos.split(':')[1]
    length = int(pos.split('-')[1]) - int(pos.split('-')[0])
    return length

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Get best match sequences", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-b", "--bed", type=str, help="bed file output from rustybam")
    parser.add_argument("-f", "--filtered", help="filtered output bed file name with sequence identity", type=str)

    args = parser.parse_args()

    bed_file = pd.read_table(args.bed)
    #bed_file = pd.read_table("NA18967_h2.bed")
    filtered_bed_file = open(args.filtered, 'w')
    bed_file['length'] = bed_file['reference_end'] - bed_file['reference_start']
    df = bed_file[['#reference_name', 'reference_start', 'reference_end', 'query_name', 'length','perID_by_matches']]
    df.columns = ['chromo','start','end','best_match','length','iden']
    df_sorted = df.sort_values(by=['chromo', 'start'])

    if df_sorted.empty:
        filtered_bed_file.close()
    else:
        print("df is not empty")
        grouped = []
        current_group = [df_sorted.iloc[0]]
        current_start = df_sorted.iloc[0]['start']
        current_end = df_sorted.iloc[0]['end']
        for i in range(1, len(df_sorted)):
            if df_sorted.iloc[i]['chromo'] == df_sorted.iloc[i-1]['chromo'] and (current_end - df_sorted.iloc[i]['start'] < (current_end - current_start) * 0.01) and (current_end - df_sorted.iloc[i]['start'] < (df_sorted.iloc[i]['end'] - df_sorted.iloc[i]['start']) * 0.01):
                grouped.append(current_group)
                current_group = [df_sorted.iloc[i]]
                current_start = min(x['start'] for x in current_group)
                current_end = max(x['end'] for x in current_group)
            elif df_sorted.iloc[i]['chromo'] == df_sorted.iloc[i-1]['chromo'] and df_sorted.iloc[i]['start'] <= current_end and df_sorted.iloc[i]['end'] >= current_start:
                current_group.append(df_sorted.iloc[i])
                current_start = min(x['start'] for x in current_group)
                current_end = max(x['end'] for x in current_group)
            else:
                grouped.append(current_group)
                current_group = [df_sorted.iloc[i]]
                current_start = min(x['start'] for x in current_group)
                current_end = max(x['end'] for x in current_group)
        grouped.append(current_group)
        ##filteration 
        sample = args.bed.split('/')[-1].split('_')[0]
        is_human = sample not in ['MFA','SSY','PPY','PAB','GGO','PPA','PTR']
        print('human' if is_human else 'primates')

        for group in grouped:
            candidates = []
            for x in group:
                ref_length = get_length(x['best_match'])
                if is_human:
                    cond = (int(x['length']) >= ref_length * 0.1) and (x['iden'] > 95)
                else:
                    cond = (ref_length * 1.2 >= int(x['length']) >= ref_length * 0.8) and (x['iden'] > 80)
                if cond:
                    x = x.copy()
                    x['ref_length'] = ref_length
                    x['ref_diff'] = abs(int(x['length']) - ref_length) / ref_length
                    candidates.append(x)

            if not candidates:
                continue

            candidates_sorted = sorted(
                candidates,
                key=lambda x: (x['iden'], int(x['length'])),
                reverse=True
            )

            if len(candidates_sorted) == 1:
                best = candidates_sorted[0]
            else:
                first = candidates_sorted[0]
                second = candidates_sorted[1]
                iden_close = abs(first['iden'] - second['iden']) < 0.1
                length_significantly_better = int(second['length']) >= int(first['length']) * 1.1
                if iden_close and length_significantly_better:
                    best = second
                else:
                    best = first
            
            #best = max(candidates, key=lambda x: int(x['length']))
            #high_iden = [c for c in candidates if c['iden'] >= 98]
            #if high_iden:
            #    best = max(high_iden, key=lambda x: int(x['length']))
            #else:
            #    best = max(candidates, key=lambda x: (x['iden'], int(x['length'])))


            outline = (
                    f"{best['chromo']}\t"
                    f"{best['start'] + 1 if best['start'] == 0 else best['start']}\t"
                    f"{best['end']}\t"
                    f"{best['best_match']}\t"
                    f"{int(best['length'])}\t"
                    f"{best['iden']}\n"
                )
            filtered_bed_file.write(outline)
