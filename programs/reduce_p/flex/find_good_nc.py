import itertools
import subprocess
from pathlib import Path

TIMEOUT_SEC = 60

# --- ここを追加 ---
Path("./good_nc").mkdir(parents=True, exist_ok=True)
# -----------------

# 0から5の数字
numbers = range(6)

# 和が奇数になるように、順序を考慮した2つの数字のペア（順列）をリスト化
valid_pairs = []
for p in itertools.permutations(numbers, 2):
    if sum(p) % 2 != 0:
        valid_pairs.append(f"{p[0]},{p[1]}")

# 1組から3組までのパターンを生成
for r in range(4, 6):
    print(f"--- {r}組のパターン ---")
    
    for combo in itertools.combinations(valid_pairs, r):
        cmd = ["./find_fg", "--p", "768", "--active", "0", "2", "4", "--outdir", "good_nc", "--save", "--nc"]
        cmd.extend(combo)
        
        print("実行中:", " ".join(cmd))
        
        # 実際にコマンドを実行する場合は以下のコメントアウトを外す
        try:
            # timeout引数を指定して実行
            subprocess.run(cmd, timeout=TIMEOUT_SEC)
        except subprocess.TimeoutExpired:
            print(f"  -> 【スキップ】処理時間が {TIMEOUT_SEC} 秒を超過しました。")
        except Exception as e:
            print(f"  -> 【エラー】予期せぬエラーが発生しました: {e}")

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