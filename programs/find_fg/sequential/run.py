import subprocess
import time

def main():
    # 探索対象のPのリスト
    p_values = [
        # 767, 766, 764, 
        763, 758, 755, 753, 752, 749, 747, 746, 745, 737, 736, 734, 731, 725, 724, 723, 722, 721, 718, 717, 716, 713, 712, 711, 707, 706, 704, 703, 699, 698, 697, 695, 694, 692, 689, 688, 687, 686, 685, 681, 679, 676, 675, 674, 671, 669, 668, 667, 664, 662, 657, 656, 655, 652, 649, 648, 640, 639, 637, 635, 634, 633, 632, 629, 628, 626, 623, 622, 621, 614, 611, 608, 605, 604, 603, 597, 596, 592, 591, 589, 586, 584, 583, 581, 579, 578, 576, 575, 573, 568, 567, 566, 565, 562, 559, 556, 554, 553, 551, 549, 548, 545, 544, 543, 542, 539, 538, 537, 536, 535, 533, 531, 527, 526, 524, 519, 517, 515, 514, 513, 511, 508, 507, 505, 502, 501, 500
    ]
    
    timeout_sec = 1200  # タイムアウト時間（秒）
    executable = "./find_fg.exe"  # コンパイル済みの実行ファイル名

    print(f"連続探索を開始する。合計 {len(p_values)} 個、各タイムアウト {timeout_sec} 秒")

    for p in p_values:
        print("\n" + "=" * 60)
        print(f"▶ P = {p} の探索を開始")
        print("=" * 60)

        start_time = time.time()
        try:
            # プロセスの実行
            # 標準出力・標準エラー出力はそのままコンソールに表示される
            result = subprocess.run(
                [executable, str(p)],
                timeout=timeout_sec
            )
            
            elapsed = time.time() - start_time
            if result.returncode == 0:
                print(f"\n[完了] P = {p} の探索が {elapsed:.1f} 秒で正常終了した。")
            else:
                print(f"\n[異常終了] P = {p} はリターンコード {result.returncode} で終了した。")

        except subprocess.TimeoutExpired:
            print(f"\n[タイムアウト] P = {p} の探索は {timeout_sec} 秒を超過したため強制終了された。次のPへ移行する。")
            
        except FileNotFoundError:
            print(f"\n[致命的エラー] 実行ファイル '{executable}' が見つからない。")
            print("事前に gcc -O3 -fopenmp find_fg_4.c -o find_fg -lm でコンパイルを行う必要がある。")
            break
            
        except KeyboardInterrupt:
            print("\n[中断] ユーザーによって全体の実行がキャンセルされた。")
            break

    print("\n全てのキューの実行が完了した。")

if __name__ == "__main__":
    main()