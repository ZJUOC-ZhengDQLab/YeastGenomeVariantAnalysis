import os
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

states_mapping = {
    '10': 'a1',
    '20': 'a2',
    '01': 'a3',
    '02': 'a4',
    '010': 'b1',
    '020': 'b2',
    '210': 'd1',
    '120': 'd3',
    '012': 'd2',
    '021': 'd4',
    '1020': 'd5',
    '2010': 'd9',
    '1010': 'd14',
    '1210': 'd55',
    '2020': 'd16',
    '2120': 'd17',
    '0201': 'd10',
    '0102': 'd6',
    '0202': 'd23',
    '0101': 'd13',
    '0121': 'd18',
    '0212': 'd25',
    '12010': 'd12',
    '12020': 'd33',
    '20210': 'd19',
    '10210': 'd67',
    '21210': 'd7',
    '21020': 'd86',
    '21010': 'd90',
    '01012': 'd15',
    '01010': 'd83',
    '02021': 'd11',
    '02012': 'd52',
    '01201': 'd21',
    '02101': 'd24',
    '02121': 'd51',
    '01212': 'd37',
    '01202': 'd29',
    '20120': 'd44',
    '102020': 'd22',
    '101020': 'd38',
    '101010': 'd42',
    '202010': 'd43',
    '201010': 'd48',
    '201020': 'd81',
    '202020': 'd68',
    '010102': 'd28',
    '010212': 'd77',
    '021201': 'd46',
    '020201': 'd54',
    '020202': 'd85',
    '020101': 'd65',
    '020121': 'd40',
    '121210': 'd47',
    '0102021': 'd26',
    '0101010': 'd84',
    '0101012': 'd53',
    '1020120': 'd20',
    '1212010': 'd35',
    '2020210': 'd36',
    '1202020': 'd39',
    '0201021': 'd30',
    '0101021': 'd82',
    '0202121': 'd59',
    '0212121': 'd76',
    '0120212': 'd49',
    '0121212': 'd66',
    '1010210': 'd79',
    '10202010': 'd27',
    '20202020': 'd62',
    '10101010': 'd64',
    '10202020': 'd70',
    '12121010': 'd71',
    '20201010': 'd31',
    '12120120': 'd41',
    '01010102': 'd32',
    '01010101': 'd72',
    '02012121': 'd69',
    '01212121': 'd78',
    '020210201': 'd61',
    '010101212': 'd50',
    '210202020': 'd73',
    '212020210': 'd80',
    '0212120202': 'd8',
    '0201212121': 'd58',
    '0212120101': 'd45',
    '0202010102': 'd74',
    '0202020102': 'd89',
    '0101010102': 'd87',
    '0101010121': 'd63',
    '01020210121': 'd60',
    '02121021210': 'd88',
    '020121212121': 'd34',
    '020101212121': 'd75',
    '021010212120212': 'd57',
    '202020201012020': 'd56',
    '0120': 'c5',
    '0210': 'c15',
    '01010': 'c7',
    '01020': 'c8',
    '02020': 'c9',
    '02010': 'c17',
    '02120': 'c4',
    '01210': 'c10',
    '021210': 'c2',
    '020210': 'c34',
    '012120': 'c16',
    '012010': 'c26',
    '010120': 'c29',
    '021020': 'c28',
    '021010': 'c35',
    '020120': 'c31',
    '012010': 'c55',
    '010120': 'c56',
    '0102010': 'c6',
    '0121020': 'c13',
    '0202020': 'c22',
    '0202010': 'c18',
    '0101010': 'c30',
    '0101020': 'c49',
    '0121210': 'c36',
    '0202120': 'c53',
    '0212010': 'c58',
    '01201010': 'c39',
    '02021210': 'c40',
    '02021010': 'c57',
    '01212010': 'c20',
    '02012010': 'c25',
    '01012010': 'c32',
    '021201020': 'c1',
    '012121010': 'c14',
    '020201010': 'c42',
    '020202020': 'c48',
    '020201210': 'c11',
    '012121210': 'c12',
    '021212020': 'c21',
    '010102020': 'c23',
    '010101020': 'c24',
    '020202010': 'c45',
    '020101010': 'c27',
    '010101010': 'c59',
    '010202010': 'c50',
    '010101210': 'c54',
    '0210101020': 'c41',
    '0212121210': 'c37',
    '0102010120': 'c51',
    '01012021010': 'c38',
    '01010101020': 'c19',
    '01010102020': 'c52',
    '020202010102': 'c46',
    '010102102121': 'c44',
    '0101021021212': 'c43',
    '0102020202010': 'c47',
    '01012101010120': 'c33'
}

event_class_mapping = {
    'a': 'CO/BIR',
    'b': 'CON',
    'd': 'CON/CO',
    'c': 'COM/CON'
}


def event_cls(data: pd.DataFrame):
    
    # cur_state: initialize -1, 0 - green, 1 - red, 2 - blue
    cur_state = -1

    history_state = ''
    events = []
    event_details = []

    last_snp = None
    last_event_snp = None
    event_points = []
    details = []

    for i in range(len(data)):
        no_hete = data.loc[i]['no_hete']
        delta = data.loc[i]['delta']
        w303 = data.loc[i]['w303']
        yjm = data.loc[i]['yjm']
        coord = data.loc[i]['coordinate']

        next_state = -1
        if not no_hete: 
            if (w303 < 0.15 and cur_state == 2):
                next_state = 2
            elif (yjm < 0.15 and cur_state == 1):
                next_state = 1
            else:
                next_state = 0
        else:
            if delta < 0.65:
                if (w303 < 0.15 and cur_state == 2):
                    next_state = 2
                elif (yjm < 0.15 and cur_state == 1):
                    next_state = 1
                else:
                    next_state = 0
            else:
                if w303 < 0.15: next_state = 2
                elif yjm < 0.15: next_state = 1
                else: next_state = 0
            
        if cur_state != -1:
            if last_event_snp is not None and coord - last_event_snp > 20000 and cur_state == 0 and history_state != '':
                history_state += str(cur_state)

                event_type = states_mapping[history_state] if history_state in states_mapping else history_state
                event_string = history_state
                details = pd.concat(details, axis=0)
                details.drop_duplicates(inplace=True)

                events.append((event_type, event_string, event_points))
                event_details.append(details)
                
                history_state = ''
                event_points = []
                details = []

            if cur_state != next_state:
                history_state += str(cur_state)
        
                cur_state = next_state

                event_points.append((last_snp, coord))

                # print(history_state, coord, cur_state, w303, yjm)

                pre_idx = i - 3 if i - 3 >= 0 else 0
                next_idx = i + 3 if i + 3 < len(data) else len(data)
                detail_ = data.loc[pre_idx:next_idx]
                detail_ = detail_[['chromosome', 'coordinate', 'w303', 'yjm']]
                detail_['sign'] = 0
                detail_.loc[detail_['coordinate'] == coord, ['sign']] = 1
                details.append(detail_)
                last_event_snp = coord
        else:
            cur_state = next_state

        last_snp = coord

    if len(event_points):
        history_state += str(cur_state)

        # a specific case
        if len(history_state) > 2 and history_state[-2:] in ['12', '21']:
            if event_points[-1][1] - event_points[-2][1] > 20000:
                # split
                history_state_1 = history_state[:-1]
                history_state_2 = history_state[-2:]

                event_type_1 = states_mapping[history_state_1] if history_state_1 in states_mapping else history_state_1
                event_type_2 = states_mapping[history_state_2] if history_state_2 in states_mapping else history_state_2
                
                event_string_1 = history_state_1
                event_string_2 = history_state_2

                event_points_1 = event_points[:-1]
                event_points_2 = event_points[-1:]

                details_1 = details[:-1]
                details_2 = details[-1:]

                details_1 = pd.concat(details_1, axis=0)
                details_2 = pd.concat(details_2, axis=0)

                details_1.drop_duplicates(inplace=True)
                details_2.drop_duplicates(inplace=True)

                events.append((event_type_1, event_string_1, event_points_1))
                events.append((event_type_2, event_string_2, event_points_2))
                event_details.append(details_1)
                event_details.append(details_2)
        else:
            event_type = states_mapping[history_state] if history_state in states_mapping else history_state
            event_string = history_state
            details = pd.concat(details, axis=0)
            details.drop_duplicates(inplace=True)
            events.append((event_type, event_string, event_points))
            event_details.append(details)

    event_details = pd.concat(event_details, axis=0) if len(event_details) else []

    return events, event_details


def output_events(result):

    all_events = []
    cur_event_id = 1

    all_event_details = []
    new_sign = []
    nan_sign = []

    for chrom, events, event_details in result:
        for item in events:
            event_type = item[0]
            event_string = item[1]
            event_class = event_class_mapping[item[0][0]] if item[0][0] in event_class_mapping else 'Unknown'

            for idx, e in enumerate(item[2]):
                all_events.append((cur_event_id, chrom, event_class, event_type, chr(ord('a') + idx), e[0], e[1]))
                new_sign.append(cur_event_id)
                if idx != 0:
                    nan_sign.append(True)
                else:
                    nan_sign.append(False)
            
            cur_event_id += 1

        all_event_details.append(event_details)

    all_events = pd.DataFrame(all_events)
    all_events.columns = ['Envent ID', 'Chromosome', 'Event type', 'LOH class', 'Transition label', 'Coordinates of LOH window (Left)', 'Coordinates of LOH window (Right)']
    all_events.loc[nan_sign, ['Envent ID', 'Chromosome', 'Event type', 'LOH class']] = np.nan
    
    all_event_details = pd.concat(all_event_details, axis=0)
    all_event_details.loc[all_event_details['sign'] == 1, ['sign']] = new_sign

    return all_events, all_event_details
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file_dir', required=True, type=str)
    parser.add_argument('--output_dir', required=True, type=str)
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)

    file_list = [os.path.join(root, f) for root, dirs, files in os.walk(args.file_dir) for f in files if f.endswith('data')]

    for file in file_list:
        filename = file.split('/')[-1]
        print(filename)

        try:
            data = pd.read_table(file, sep=" ")
            data.columns = ['chromosome', 'coordinate', 'w303', 'yjm']
        except Exception as e:
            print(f"Error reading file {filename}: {e}")
            continue

        if len(data) == 0 or 'chromosome' not in data.columns:
            print(f"Skipping file {filename} due to invalid format or empty data.")
            continue

        chrom_data = []
        for chrom in data['chromosome'].unique():
            sub_data = data[data['chromosome'] == chrom].reset_index(drop=True)
            sub_data['w303'] = np.clip(sub_data['w303'], 0.0, 1.0)
            sub_data['yjm'] = np.clip(sub_data['yjm'], 0.0, 1.0)

            sub_data['hete_1'] = (0.5 - np.abs(sub_data['w303'] - 0.5)) / 0.5
            sub_data['hete_2'] = (0.5 - np.abs(sub_data['yjm'] - 0.5)) / 0.5

            sub_data['no_hete_1'] = 1 - sub_data['hete_1']
            sub_data['no_hete_2'] = 1 - sub_data['hete_2']
            sub_data['no_hete'] = (sub_data['no_hete_1'] >= 0.90) | (sub_data['no_hete_2'] >= 0.90)

            sub_data['delta'] = np.abs(sub_data['w303'] - sub_data['yjm']) / 1.0

            events, event_details = event_cls(sub_data)
            if len(events):
                chrom_data.append((chrom, events, event_details))

        if len(chrom_data):
            all_events, all_event_details = output_events(chrom_data)
            with pd.ExcelWriter(os.path.join(args.output_dir, f'{filename}.xlsx')) as writer:
                all_events.to_excel(writer, sheet_name='Events', index=False)
                all_event_details.to_excel(writer, sheet_name='Details', index=False)


if __name__ == '__main__':
    main()