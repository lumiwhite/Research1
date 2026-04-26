from pathlib import Path

# 対象となるディレクトリの設定
target_dir = Path("results")

# 検索する文字列
search_string = "failures=0"

# ディレクトリが存在するか確認
if not target_dir.exists() or not target_dir.is_dir():
    print(f"エラー: '{target_dir}' ディレクトリが見つからない。")
else:
    # ディレクトリ直下のアイテムを順番に確認
    for filepath in target_dir.iterdir():
        # 対象がファイルである場合のみ処理（ディレクトリ等はスキップ）
        if filepath.is_file():
            try:
                # ファイルを開いて文字列を検索する
                # 文字コードのエラーで停止しないよう errors="ignore" を付与
                with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                    content = f.read()
                    if search_string in content:
                        print(filepath.name)
            except Exception as e:
                # 権限エラー等で読み込めなかった場合の例外処理
                print(f"警告: {filepath.name} の読み込み中にエラーが発生した ({e})")