#!/bin/bash

# 定义参数
parameter="1000000"

# 获取输入文件目录和输出文件目录
input_dir="./bin_data"
output_dir="./test_result"

# 确保输出目录存在
mkdir -p "$output_dir"

# 遍历输入目录中的所有二进制文件
for input_file in "$input_dir"/*.bin; do
  # 获取输入文件名（不带路径和扩展名）
  base_filename=$(basename "$input_file" .bin)

  # 定义生成文件的目标路径和名称
  output_file="$output_dir/${base_filename}_nist_output.txt"

  # 运行 Expect 脚本来执行 NIST 测试程序
  expect ./nist_interact.exp "$input_file" "$parameter"

  # 确保 Expect 脚本执行结束
  wait

  # 假设程序生成的输出文件是 generated_output_file.txt
  generated_file="./sts-2.1.2_optimized/experiments/AlgorithmTesting/finalAnalysisReport.txt"

  # 检查文件是否生成
  if [ -f "$generated_file" ]; then
    # 将生成的文件移动并重命名
    mv "$generated_file" "$output_file"
    echo "Moved and renamed output to: $output_file"
  else
    echo "Error: Generated file not found for $input_file!"
  fi
done
