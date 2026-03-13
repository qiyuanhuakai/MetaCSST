#!/usr/bin/env python3
"""
MetaCSST - Numba optimized version with pyfastx lazy loading
该文件中的search_vr_core和call_vr函数的if语句嵌套过于严重，导致代码可读性下降。我期待能在未来修复它，但是由于numba的特性，这也需要设计复杂的参数传入传出
"""

import sys
import logging
from typing import Tuple, Dict, List

import numpy as np
from numpy.typing import NDArray
from numba import njit
import pyfastx

# 配置日志格式
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stderr
)
logger = logging.getLogger(__name__)

# 全局常量
NA_RATIO: float = 0.5
MISS: int = 3
LEN_MIN: int = 30
IDENTITY: float = 0.6


# 碱基编码: A=0, T=1, C=2, G=3, other=4
@njit(cache=True)
def encode_base(c: int) -> int:
    """
    将单个碱基的ASCII码转换为编码值。
    
    Args:
        c: 碱基字符的ASCII码
        
    Returns:
        编码值: A=0, T=1, C=2, G=3, other=4
    """
    if c == 65 or c == 97:  # A, a
        return 0
    elif c == 84 or c == 116:  # T, t
        return 1
    elif c == 67 or c == 99:  # C, c
        return 2
    elif c == 71 or c == 103:  # G, g
        return 3
    return 4


@njit(cache=True)
def encode_seq(seq_bytes: NDArray[np.uint8]) -> NDArray[np.int8]:
    """
    将序列的字节数组转换为编码数组。
    
    Args:
        seq_bytes: 序列的uint8字节数组
        
    Returns:
        编码后的int8数组
    """
    n = len(seq_bytes)
    arr = np.empty(n, dtype=np.int8)
    for i in range(n):
        arr[i] = encode_base(seq_bytes[i])
    return arr


@njit(cache=True)
def reverse_complement(arr: NDArray[np.int8]) -> NDArray[np.int8]:
    """
    计算序列的反向互补序列。
    
    反向互补规则: A<->T(0<->1), C<->G(2<->3)
    
    Args:
        arr: 编码后的序列数组
        
    Returns:
        反向互补后的序列数组
    """
    n = len(arr)
    result = np.empty(n, dtype=np.int8)
    for i in range(n):
        c = arr[n - 1 - i]
        if c == 0:
            result[i] = 1
        elif c == 1:
            result[i] = 0
        elif c == 2:
            result[i] = 3
        elif c == 3:
            result[i] = 2
        else:
            result[i] = c
    return result


@njit(cache=True)
def num_NA(arr: NDArray[np.int8]) -> int:
    """
    计算序列中非A碱基的数量。
    
    Args:
        arr: 编码后的序列数组
        
    Returns:
        非A碱基的计数
    """
    count = 0
    for i in range(len(arr)):
        if arr[i] != 0:
            count += 1
    return count


@njit(cache=True)
def re_check_numba(
    TR_arr: NDArray[np.int8],
    VR_arr: NDArray[np.int8],
    mut_A: int,
    mut_NA: int,
    identity_threshold: float
) -> bool:
    """
    验证TR-VR配对是否符合条件。
    
    Args:
        TR_arr: TR序列的编码数组
        VR_arr: VR序列的编码数组
        mut_A: 预期的A碱基突变数
        mut_NA: 预期的非A碱基错误数
        identity_threshold: 相似度阈值
        
    Returns:
        配对是否有效
    """
    n = len(TR_arr)
    mut = 0
    err = 0
    same = 0
    for i in range(n):
        c1, c2 = TR_arr[i], VR_arr[i]
        if c1 != c2:
            if c1 == 0:  # A
                mut += 1
            else:
                err += 1
        else:
            same += 1

    # 模拟Perl的sprintf("%0.2f", ...)舍入
    identity = round(same / n * 100) / 100.0

    if mut == mut_A and err == mut_NA and identity >= identity_threshold:
        return True
    return False


@njit(cache=True)
def search_VR_core(
    TR_arr: NDArray[np.int8],
    start_param: int,
    end_param: int,
    seq_arr: NDArray[np.int8],
    seq2_arr: NDArray[np.int8],
    seq1_arr: NDArray[np.int8],
    miss: int,
    is_positive_strand: bool,
    identity: float,
    len_min: int,
) -> Tuple[NDArray[np.int64], int, int, int]:
    """
    核心搜索函数，在序列中搜索VR区域。
    
    Args:
        TR_arr: TR序列的编码数组
        start_param: 搜索起始位置
        end_param: 搜索结束位置
        seq_arr: 原始序列的编码数组
        seq2_arr: 正链或反向互补序列的编码数组
        seq1_arr: 反向互补序列的编码数组
        miss: 允许的最大错配数
        is_positive_strand: 是否为正链
        identity: 相似度阈值
        len_min: 最小长度要求
        
    Returns:
        元组包含:
        - 结果数组 (每行: [a, b, c, d, mut, error, is_negative_vr])
        - 结果数量
        - total_start
        - total_end
    """
    tr_len = len(TR_arr)
    seq_len = len(seq_arr)

    # 预分配结果数组（最多1000个结果）
    max_results = 1000
    results = np.zeros((max_results, 7), dtype=np.int64)
    result_count = 0

    total_start = start_param
    total_end = end_param

    if is_positive_strand:
        start = start_param
        end = end_param
    else:
        start = seq_len - 1 - end_param
        end = seq_len - 1 - start_param

    # Search positive strand VR
    for i in range(seq_len - tr_len):
        if i + tr_len - 1 < start_param or i > end_param:
            error = 0
            mut = 0
            skip = False

            for j in range(tr_len):
                c1 = TR_arr[j]
                c2 = seq_arr[i + j]
                if c1 == 0 and c1 != c2:  # A
                    mut += 1
                elif c1 != 0 and c1 != c2:
                    error += 1
                    if error > miss:
                        skip = True
                        break

            if not skip:
                # 左扩展
                left = 1
                while i - left >= 0 and start - left >= 0:
                    c1 = seq2_arr[start - left]
                    c2 = seq_arr[i - left]
                    if c1 == 0 and c1 != c2:
                        mut += 1
                    elif c1 != 0 and c1 != c2:
                        break
                    left += 1

                # 右扩展
                right = 1
                while i + tr_len - 1 + right < seq_len and end + right < seq_len:
                    c1 = seq2_arr[end + right]
                    c2 = seq_arr[i + tr_len - 1 + right]
                    if c1 == 0 and c1 != c2:
                        mut += 1
                    elif c1 != 0 and c1 != c2:
                        break
                    right += 1

                left -= 1
                right -= 1

                a = start - left
                b = end + right
                c = i - left
                d = i + tr_len - 1 + right

                tr_new_len = b - a + 1
                vr_new_len = d - c + 1

                if tr_new_len == vr_new_len and tr_new_len >= len_min:
                    TR_new = seq2_arr[a : b + 1]
                    VR_new = seq_arr[c : d + 1]

                    if re_check_numba(TR_new, VR_new, mut, error, identity):
                        # 坐标转换
                        a_out, b_out = a, b
                        if not is_positive_strand:
                            a_out = seq_len - 1 - b
                            b_out = seq_len - 1 - a

                        if result_count < max_results:
                            results[result_count, 0] = a_out
                            results[result_count, 1] = b_out
                            results[result_count, 2] = c
                            results[result_count, 3] = d
                            results[result_count, 4] = mut
                            results[result_count, 5] = error
                            results[result_count, 6] = 0  # positive VR
                            result_count += 1

                        if a_out < total_start:
                            total_start = a_out
                        if c < total_start:
                            total_start = c
                        if b_out > total_end:
                            total_end = b_out
                        if d > total_end:
                            total_end = d

    # Search negative strand VR
    for i in range(seq_len - tr_len):
        if i + tr_len - 1 < seq_len - 1 - end_param or i > seq_len - 1 - start_param:
            error = 0
            mut = 0
            skip = False

            for j in range(tr_len):
                c1 = TR_arr[j]
                c2 = seq1_arr[i + j]
                if c1 == 0 and c1 != c2:  # A
                    mut += 1
                elif c1 != 0 and c1 != c2:
                    error += 1
                    if error > miss:
                        skip = True
                        break

            if not skip:
                left = 1
                while i - left >= 0 and start - left >= 0:
                    c1 = seq2_arr[start - left]
                    c2 = seq1_arr[i - left]
                    if c1 == 0 and c1 != c2:
                        mut += 1
                    elif c1 != 0 and c1 != c2:
                        break
                    left += 1

                right = 1
                while i + tr_len - 1 + right < seq_len and end + right < seq_len:
                    c1 = seq2_arr[end + right]
                    c2 = seq1_arr[i + tr_len - 1 + right]
                    if c1 == 0 and c1 != c2:
                        mut += 1
                    elif c1 != 0 and c1 != c2:
                        break
                    right += 1

                left -= 1
                right -= 1

                a = start - left
                b = end + right
                c = i - left
                d = i + tr_len - 1 + right

                tr_new_len = b - a + 1
                vr_new_len = d - c + 1

                if tr_new_len == vr_new_len and tr_new_len >= len_min:
                    TR_new = seq2_arr[a : b + 1]
                    VR_new = seq1_arr[c : d + 1]

                    if re_check_numba(TR_new, VR_new, mut, error, identity):
                        # VR坐标转换
                        c_out = seq_len - 1 - d
                        d_out = seq_len - 1 - c

                        # TR坐标转换
                        a_out, b_out = a, b
                        if not is_positive_strand:
                            a_out = seq_len - 1 - b
                            b_out = seq_len - 1 - a

                        if result_count < max_results:
                            results[result_count, 0] = a_out
                            results[result_count, 1] = b_out
                            results[result_count, 2] = c_out
                            results[result_count, 3] = d_out
                            results[result_count, 4] = mut
                            results[result_count, 5] = error
                            results[result_count, 6] = 1  # negative VR
                            result_count += 1

                        if a_out < total_start:
                            total_start = a_out
                        if c_out < total_start:
                            total_start = c_out
                        if b_out > total_end:
                            total_end = b_out
                        if d_out > total_end:
                            total_end = d_out

    return results[:result_count], result_count, total_start, total_end


def decode_seq(arr: NDArray[np.int8]) -> str:
    """
    将编码数组转回碱基字符串。
    
    Args:
        arr: 编码后的序列数组
        
    Returns:
        碱基序列字符串
    """
    bases = ["A", "T", "C", "G", "N"]
    return "".join(bases[c] if c < 5 else "N" for c in arr)


def my_reverse(seq: str) -> str:
    """
    计算序列的反向互补（Python字符串版本）。
    
    Args:
        seq: DNA序列字符串
        
    Returns:
        反向互补序列字符串
    """
    complement: Dict[str, str] = {
        "A": "T", "a": "T",
        "T": "A", "t": "A",
        "C": "G", "c": "G",
        "G": "C", "g": "C",
    }
    return "".join(complement.get(c, c) for c in reversed(seq))


def search_VR(
    TR: str,
    start_param: int,
    end_param: int,
    seq: str,
    miss: int,
    string: str,
    seq_id: str,
    identity: float
) -> Tuple[str, int, int, int]:
    """
    封装函数，调用numba优化的核心搜索。
    
    Args:
        TR: TR序列字符串
        start_param: 起始位置
        end_param: 结束位置
        seq: 完整基因组序列
        miss: 允许的错配数
        string: 链方向 ("+" 或 "-")
        seq_id: 序列ID
        identity: 相似度阈值
        
    Returns:
        元组包含:
        - 输出字符串
        - 找到的VR数量
        - total_start
        - total_end
    """

    seq_len = len(seq)

    # 编码序列
    TR_arr = encode_seq(np.frombuffer(TR.encode(), dtype=np.uint8))
    seq_arr = encode_seq(np.frombuffer(seq.encode(), dtype=np.uint8))
    seq1_arr = reverse_complement(seq_arr)  # 反向互补

    is_positive = string == "+"
    if is_positive:
        seq2_arr = seq_arr.copy()
    else:
        seq2_arr = reverse_complement(seq_arr)

    # 调用核心函数
    results, result_count, total_start, total_end = search_VR_core(
        TR_arr,
        start_param,
        end_param,
        seq_arr,
        seq2_arr,
        seq1_arr,
        miss,
        is_positive,
        identity,
        LEN_MIN,
    )

    # 生成输出字符串
    tmp = ""
    if result_count > 0:
        # 预计算seq2用于提取序列
        seq2 = seq if is_positive else my_reverse(seq)
        seq1 = my_reverse(seq)

        for idx in range(result_count):
            a, b, c, d = (
                int(results[idx, 0]),
                int(results[idx, 1]),
                int(results[idx, 2]),
                int(results[idx, 3]),
            )
            mut, error = int(results[idx, 4]), int(results[idx, 5])
            is_neg_vr = results[idx, 6] == 1

            # 重新计算内部坐标来提取序列
            if is_positive:
                a_internal, b_internal = a, b
            else:
                a_internal = seq_len - 1 - b
                b_internal = seq_len - 1 - a

            TR_new = seq2[a_internal : b_internal + 1]

            if is_neg_vr:
                c_internal = seq_len - 1 - d
                d_internal = seq_len - 1 - c
                VR_new = seq1[c_internal : d_internal + 1]
                vr_strand = "-"
            else:
                VR_new = seq[c : d + 1]
                vr_strand = "+"

            tmp += f"{seq_id}\tTR\t{string}\t{start_param}\t{end_param}\t{a}\t{b}\t{mut}\t{error}\t{TR_new}\n"
            tmp += f"{seq_id}\tVR\t{vr_strand}\t*\t*\t{c}\t{d}\t{mut}\t{error}\t{VR_new}\n"

    return tmp, result_count, total_start, total_end


def call_vr(predict_file: str, fasta_index: pyfastx.Fasta) -> str:
    """
    处理预测文件，调用VR搜索。
    
    Args:
        predict_file: 输入的GTF预测文件路径
        fasta_index: pyfastx.Fasta索引对象
        
    Returns:
        处理后的输出字符串
    """
    output = ""
    miss_half = MISS // 2

    sub_seq: List[str] = []
    sub_start: List[int] = []
    sub_end: List[int] = []
    sub_string: List[str] = []
    type_list: List[int] = []

    seq_id = ""
    seq = ""
    start = end = pair = 0
    tmp = ""

    # 缓存已读取的序列，避免重复读取同一序列
    seq_cache: Dict[str, str] = {}

    with open(predict_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue

            fields = line.split()
            if "TR" in line:
                sub_seq.append(fields[6])
                sub_start.append(int(fields[4]))
                sub_end.append(int(fields[5]))
                sub_string.append(fields[3])
                type_list.append(1)
            elif "VR" in line:
                sub_seq.append(fields[6])
                sub_start.append(int(fields[4]))
                sub_end.append(int(fields[5]))
                sub_string.append(fields[3])
                type_list.append(2)
            elif "RT" in line:
                sub_seq.append(fields[6])
                sub_start.append(int(fields[4]))
                sub_end.append(int(fields[5]))
                sub_string.append(fields[3])
                type_list.append(3)
            elif "DGR" in line:
                start = int(fields[4])
                end = int(fields[5])
                seq_id = fields[0]

                # 按需从pyfastx获取序列（带缓存）
                if seq_id in seq_cache:
                    seq = seq_cache[seq_id]
                elif seq_id in fasta_index:
                    seq = str(fasta_index[seq_id].seq)
                    seq_cache[seq_id] = seq  # 缓存以备后续使用
                else:
                    logger.warning(f"Sequence '{seq_id}' not found in FASTA index, skipping.")
                    seq = ""

                if seq:
                    for i in range(len(sub_seq)):
                        if type_list[i] == 1:
                            start_new = -1
                            end_new = -1
                            tmp_sub = ""
                            index = 0

                            NA = int(len(sub_seq[i]) * NA_RATIO)
                            NA_count = sum(1 for c in sub_seq[i] if c != "A")

                            if NA_count >= NA:
                                quarter = (sub_end[i] - sub_start[i]) // 4 + sub_start[i]
                                middle = ((sub_end[i] - sub_start[i]) // 4) * 2 + sub_start[i]
                                quarter3 = ((sub_end[i] - sub_start[i]) // 4) * 3 + sub_start[i]

                                if sub_string[i] == "+":
                                    seq_left = seq[sub_start[i] : middle + 1]
                                    seq_middle = seq[quarter : quarter3 + 1]
                                    seq_right = seq[middle : sub_end[i] + 1]

                                    tmp_sub, index, start_new, end_new = search_VR(
                                        sub_seq[i], sub_start[i], sub_end[i],
                                        seq, MISS, sub_string[i], seq_id, IDENTITY,
                                    )

                                    if index == 0:
                                        tmp_sub, index, start_new, end_new = search_VR(
                                            seq_right, middle, sub_end[i],
                                            seq, miss_half, sub_string[i], seq_id, IDENTITY,
                                        )
                                        if index == 0:
                                            tmp_sub, index, start_new, end_new = search_VR(
                                                seq_middle, quarter, quarter3,
                                                seq, miss_half, sub_string[i], seq_id, IDENTITY,
                                            )
                                            if index == 0:
                                                tmp_sub, index, start_new, end_new = search_VR(
                                                    seq_left, sub_start[i], middle,
                                                    seq, miss_half, sub_string[i], seq_id, IDENTITY,
                                                )
                                else:
                                    seq_left = my_reverse(seq[middle : sub_end[i] + 1])
                                    seq_middle = my_reverse(seq[quarter : quarter3 + 1])
                                    seq_right = my_reverse(seq[sub_start[i] : middle + 1])

                                    tmp_sub, index, start_new, end_new = search_VR(
                                        sub_seq[i], sub_start[i], sub_end[i],
                                        seq, MISS, sub_string[i], seq_id, IDENTITY,
                                    )

                                    if index == 0:
                                        tmp_sub, index, start_new, end_new = search_VR(
                                            seq_right, sub_start[i], middle,
                                            seq, miss_half, sub_string[i], seq_id, IDENTITY,
                                        )
                                        if index == 0:
                                            tmp_sub, index, start_new, end_new = search_VR(
                                                seq_middle, quarter, quarter3,
                                                seq, miss_half, sub_string[i], seq_id, IDENTITY,
                                            )
                                            if index == 0:
                                                tmp_sub, index, start_new, end_new = search_VR(
                                                    seq_left, middle, sub_end[i],
                                                    seq, miss_half, sub_string[i], seq_id, IDENTITY,
                                                )

                            tmp += tmp_sub
                            pair += index

                            if start_new != -1 and start_new < start:
                                start = start_new
                            if end_new > end:
                                end = end_new

                    if pair != 0:
                        output += tmp
                        for i in range(len(sub_seq)):
                            if type_list[i] == 3:
                                output += f"{seq_id}\tRT\t{sub_string[i]}\t{sub_start[i]}\t{sub_end[i]}\t{sub_seq[i]}\n"

                # Reset
                sub_seq = []
                sub_start = []
                sub_end = []
                sub_string = []
                type_list = []
                seq_id = ""
                seq = ""
                start = end = pair = 0
                tmp = ""

    return output


def remove_repeat(input_text: str) -> str:
    """
    去除重复的TR-VR配对。
    
    Args:
        input_text: 包含TR/VR/RT记录的输入文本
        
    Returns:
        去重后的输出文本
    """
    output = ""
    seen: Dict[str, int] = {}
    TR_line = ""

    for line in input_text.split("\n"):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 2:
            continue

        if fields[1] == "TR":
            TR_line = line
        elif fields[1] == "VR":
            TR_fields = TR_line.split()
            VR_fields = fields
            key = (
                TR_fields[0] + TR_fields[1] + TR_fields[2] + TR_fields[5] + TR_fields[6]
                + VR_fields[1] + VR_fields[2] + VR_fields[5] + VR_fields[6]
            )
            if key not in seen:
                output += TR_line + "\n" + line + "\n"
                seen[key] = 1
        else:
            output += line + "\n"

    return output


def main() -> None:
    """主函数：解析参数并执行VR搜索流程。"""
    if len(sys.argv) != 4:
        logger.error(
            f"Usage: {sys.argv[0]} <raw.gtf> <genome.fasta> <output.gtf>"
        )
        sys.exit(1)

    input_gtf: str = sys.argv[1]
    input_fasta: str = sys.argv[2]
    output_gtf: str = sys.argv[3]

    # 使用pyfastx建立索引（首次会创建.fxi索引文件，后续直接使用）
    logger.info("Building/loading pyfastx index...")
    fasta_index: pyfastx.Fasta = pyfastx.Fasta(input_fasta, build_index=True)
    logger.info(f"Indexed {len(fasta_index)} sequences")

    # 预热JIT编译
    logger.info("Warming up JIT compiler...")
    dummy_seq = "ATCGATCGATCG" * 10
    dummy_arr = encode_seq(np.frombuffer(dummy_seq.encode(), dtype=np.uint8))
    _ = reverse_complement(dummy_arr)
    logger.info("JIT warm-up complete")

    # 传递pyfastx索引对象，按需读取序列
    logger.info(f"Processing input file: {input_gtf}")
    vr_output: str = call_vr(input_gtf, fasta_index)
    
    logger.info("Removing duplicate entries...")
    final_output: str = remove_repeat(vr_output)

    logger.info(f"Writing output to: {output_gtf}")
    with open(output_gtf, "w") as f:
        f.write(final_output)

    logger.info("Done.")


if __name__ == "__main__":
    main()
