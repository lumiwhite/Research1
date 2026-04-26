import itertools
import subprocess
from pathlib import Path
# 1. ディレクトリの設定
input_dir = Path("./good_nc")
output_dir = Path("./results")

# 出力先ディレクトリを作成（既に存在する場合は何もしない）
output_dir.mkdir(parents=True, exist_ok=True)

# 2. ディレクトリの存在確認とループ処理
if not input_dir.exists() or not input_dir.is_dir():
    print(f"エラー: {input_dir} ディレクトリが見つからない。")
else:
    # input_dirの中身を一つずつ確認
    for filepath in input_dir.iterdir():
        # サブディレクトリ等をスキップし、ファイルのみを対象とする
        if filepath.is_file():
            filename = filepath.name
            output_filepath = output_dir / f"out_{filename}.txt"
            
            print(f"実行中: {filename} ...")
            
            # 実行するコマンドをリスト形式で定義
            cmd = [
                "./jointbp_llr",
                "--simulate",
                "--p", "0.04000000",
                "--trials", "100",
                "--max-iter", "1000",
                "--flip-hist", "5",
                "--damping", "0.000000",
                "--report-every", "20",
                "--no-pp",
                "--params", str(filepath)
            ]
            
            # コマンドを実行し、標準出力をファイルに書き込む
            with open(output_filepath, "w", encoding="utf-8") as f:
                # stdout=f とすることで、ターミナルへの出力をファイルへリダイレクトする
                subprocess.run(cmd, stdout=f, text=True)
                
    print("すべての処理が完了した。")