import sys

def fix_liftoff_duplicates(gff_file, output_file):
    id_tracker = {}
    parent_map = {}

    with open(gff_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            
            parts = line.strip().split('\t')
            if len(parts) != 9:
                outfile.write(line)
                continue
            
            attr_str = parts[8]
            attr_dict = {}
            for attr in attr_str.strip().split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            
            # 1. Parent 속성 업데이트 (자식 Feature인 경우)
            # 부모가 이전 단계에서 이름이 바뀌었다면, 자식의 Parent 정보도 바꿔야 함
            if 'Parent' in attr_dict:
                original_parents = attr_dict['Parent'].split(',')
                new_parents = []
                for p in original_parents:
                    # 가장 최근에 매핑된 부모 ID를 사용 (같은 ID의 유전자가 여러 개일 경우 뒤쪽 카피를 따라감)
                    if p in parent_map:
                        new_parents.append(parent_map[p])
                    else:
                        new_parents.append(p)
                attr_dict['Parent'] = ','.join(new_parents)

            # 2. ID 속성 업데이트 (중복 체크 및 재할당)
            if 'ID' in attr_dict:
                original_id = attr_dict['ID']
                
                # 이 ID가 처음 등장했는지, 또 등장했는지 확인
                if original_id not in id_tracker:
                    id_tracker[original_id] = 1
                    # 첫 등장이면 매핑 테이블에 원본 그대로 등록
                    parent_map[original_id] = original_id
                else:
                    # 이미 등장한 ID라면 카운트 증가시키고 suffix 부여
                    # 단, 자식 feature(예: exon, CDS)는 부모가 바뀌지 않았다면 굳이 ID를 안 바꿀 수도 있지만
                    # Liftoff는 보통 통째로 카피하므로 안전하게 모두 바꿈
                    
                    # 여기서 중요한 점: 단순 하위 feature의 연속(exon 1, exon 2...)인지
                    # 아니면 진짜 새로운 copy 유전자의 시작인지 구분해야 함.
                    # Liftoff 결과에서 'gene' 타입이 중복되면 확실히 새로운 copy임.
                    feature_type = parts[2]
                    
                    # 유전자 레벨이거나, 부모가 바뀌어서 자식도 바뀌어야 하는 경우
                    if feature_type in ['gene', 'mRNA', 'transcript']:
                        id_tracker[original_id] += 1
                        new_id = f"{original_id}_copy{id_tracker[original_id]}"
                        attr_dict['ID'] = new_id
                        parent_map[original_id] = new_id # 자식들이 이 새 ID를 참조하도록 맵 업데이트
                    else:
                        # Exon/CDS 등은 ID가 겹쳐도 보통 Parent로 구분되므로 
                        # 부모가 바뀌었으면 따라서 ID를 바꿔줌
                        if 'Parent' in attr_dict:
                             # 부모가 매핑된 상태라면 나도 바꾼다.
                             original_parent = attr_dict['Parent'].split(',')[0] # 원본 Parent 값 필요하지만 위에서 덮어씌워짐.
                             # 로직 단순화: 그냥 ID가 중복되면 무조건 copy를 붙임 (안전함)
                             id_tracker[original_id] += 1
                             new_id = f"{original_id}_copy{id_tracker[original_id]}"
                             attr_dict['ID'] = new_id
                             parent_map[original_id] = new_id

            # 3. 파일 쓰기
            new_attr_parts = [f"{k}={v}" for k, v in attr_dict.items()]
            parts[8] = ";".join(new_attr_parts)
            outfile.write('\t'.join(parts) + '\n')

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python fix_liftoff.py <input.gff3> <output.gff3>")
    else:
        fix_liftoff_duplicates(sys.argv[1], sys.argv[2])