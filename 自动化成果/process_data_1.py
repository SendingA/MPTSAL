import glob
import os


def process_line(line: str) -> str:
    """处理单行文本：移除"00> "前缀，清理空行和多余空格"""
    if line.strip().startswith("00> "):
        processed = line.strip()[4:]  # 跳过"00> "前缀（含空格，共4字符）
        return processed
    else:
        if line.strip().strip("00>"):
            processed = line.strip()[3:]
            return processed
    return line.strip()


def is_block_complete(block_content: list[str]) -> bool:
    """判定单个数据包是否完整（保留原校验逻辑）"""
    # 条件1：包含"Printing STS0 CIR"或"Printing STS1 CIR"标识
    has_cir_flag = any(
        line.startswith("Printing STS0 CIR") or line.startswith("Printing STS1 CIR")
        for line in block_content
    )

    # 条件2：包含至少1组数值数据（匹配"x,y,"格式）
    has_numeric_data = True

    # 条件3：以"_________________________________"分隔线结尾
    ends_with_separator = True
    print(f"条件是:has_cir_flag:{has_cir_flag},numeric_data:{has_numeric_data},ends_with_separator:{ends_with_separator}")
    return has_cir_flag and has_numeric_data and ends_with_separator


def parse_and_filter_complete_blocks(file_path: str) -> list[dict]:
    """解析文件并过滤不完整数据包"""
    seq_blocks = []
    current_block = None

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            processed_line = process_line(line)
            if not processed_line:
                continue

            # 识别Sequence Number行，开启新数据块
            if processed_line.startswith("Sequence Number: "):
                seq_str = processed_line.split(": ")[1]
                seq = int(seq_str)

                # 校验并保存上一个数据块
                if current_block is not None:
                    if is_block_complete(current_block["content"]):
                        seq_blocks.append(current_block)
                        print(f"  保留完整数据包：Sequence Number = {current_block['seq']}")
                    else:
                        print(f"  删除不完整数据包：Sequence Number = {current_block['seq']}")

                # 初始化新数据块
                current_block = {
                    "seq": seq,
                    "content": [processed_line]
                }

            # 非序号行，加入当前数据块
            elif current_block is not None:
                current_block["content"].append(processed_line)

    # 处理最后一个数据块
    if current_block is not None:
        if is_block_complete(current_block["content"]):
            seq_blocks.append(current_block)
            print(f"  保留完整数据包：Sequence Number = {current_block['seq']}")
        else:
            print(f"  删除不完整数据包：Sequence Number = {current_block['seq']}")

    # 无完整数据包时抛出异常
    if not seq_blocks:
        raise ValueError("文件中所有数据包均不完整，无可用数据")
    if len(seq_blocks) < 2:
        raise ValueError(f"完整数据包仅{len(seq_blocks)}个，不足2个，无法处理")

    return seq_blocks


def filter_by_head_tail(complete_blocks: list[dict]) -> list[dict]:
    """按规则删除开头和结尾的数据包"""
    blocks = complete_blocks.copy()
    n = len(blocks)
    print(f"  初始完整数据包数量：{n}，序号：{[b['seq'] for b in blocks]}")

    # 处理开头删除
    head_delete_count = 0
    for i in range(n - 1):
        current_seq = blocks[i]["seq"]
        next_seq = blocks[i + 1]["seq"]
        if current_seq + 3 == next_seq:
            blocks = blocks[i:]
            print(f"  开头删除{head_delete_count}个包，剩余序号：{[b['seq'] for b in blocks]}")
            break
        head_delete_count += 1
    else:
        raise ValueError("开头未找到满足「当前序号+3=下一序号」的数据包，无法停止删除")

    if len(blocks) < 2:
        raise ValueError(f"开头处理后仅剩余{len(blocks)}个包，不足2个")

    # 处理结尾删除
    tail_delete_count = 1
    blocks = blocks[:-1]
    print(f"  结尾先删除1个包（最后一个），剩余序号：{[b['seq'] for b in blocks]}")

    if len(blocks) < 2:
        raise ValueError("结尾删除最后一个包后，剩余不足2个包，无法继续处理")

    while len(blocks) >= 2:
        last_idx = len(blocks) - 1
        current_seq = blocks[last_idx]["seq"]
        prev_seq = blocks[last_idx - 1]["seq"]
        if current_seq - 3 == prev_seq:
            print(f"  结尾共删除{tail_delete_count}个包，剩余序号：{[b['seq'] for b in blocks]}")
            break
        blocks = blocks[:-1]
        tail_delete_count += 1
    else:
        raise ValueError("结尾未找到满足「当前序号-3=上一序号」的数据包，无法停止删除")

    if len(blocks) < 2:
        raise ValueError(f"结尾处理后仅剩余{len(blocks)}个包，不足2个")

    return blocks


def generate_output(valid_blocks: list[dict]) -> str:
    """生成最终输出内容"""
    output_lines = []
    for block in valid_blocks:
        output_lines.extend(block["content"])
        output_lines.append("_________________________________")

    if output_lines and output_lines[-1] == "_________________________________":
        output_lines.pop()

    return "\n".join(output_lines)


def batch_process(input_pattern: str, output_dir: str = "processed_files") -> None:
    """批量处理文件，输出为.txt格式"""
    os.makedirs(output_dir, exist_ok=True)
    input_files = glob.glob(input_pattern)

    if not input_files:
        print(f"未找到匹配的输入文件（模式：{input_pattern}）")
        return

    for file_path in input_files:
        file_name = os.path.basename(file_path)
        print(f"\n========================================\n")
        print(f"开始处理文件：{file_name}")
        try:
            complete_blocks = parse_and_filter_complete_blocks(file_path)
            valid_blocks = filter_by_head_tail(complete_blocks)
            print(f"  最终保留的数据包数量：{len(valid_blocks)}")

            output_content = generate_output(valid_blocks)
            # 修改输出文件后缀为.txt
            base_name = os.path.splitext(file_name)[0]
            output_file = os.path.join(output_dir, f"processed_{base_name}.txt")
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(output_content)

            print(f"\n  处理成功！输出文件：{output_file}")
        except Exception as e:
            print(f"\n  处理失败：{str(e)}")


if __name__ == "__main__":
    INPUT_FILE_PATTERN = "test-30_*.log"
    OUTPUT_DIRECTORY = "processed_rtt_files"

    print(f"批量处理开始，输入模式：{INPUT_FILE_PATTERN}，输出目录：{OUTPUT_DIRECTORY}")
    print(f"处理规则：\n  - 完全删除含00>的行\n  - 开头删除至「当前序号+3=下一序号」\n  - 结尾先删最后1个，再删除至「当前序号-3=上一序号」\n  - 输出格式为.txt\n")
    batch_process(INPUT_FILE_PATTERN, OUTPUT_DIRECTORY)
    print(f"\n========================================\n")
    print("批量处理结束！")