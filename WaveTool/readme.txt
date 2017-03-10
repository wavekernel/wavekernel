実行コマンド
python Main_AutoWave.py <elses-generate-cubefile> <JSON> <XYZ> <BASIS>
例)
python Main_AutoWave.py elses-generate-cubefile out.json position.xyz output_basis_information.txt -alpha 2.0 -s 10 --parallel --position
mpi4pyをインストールしていてmpiで並列するとき
mpiexec -n 2 python Main_AutoWave.py elses-generate-cubefile out.json position.xyz output_basis_information.txt -alpha 2.0 -s 10 --parallel --position

オプション
-h : help
-s STRIDE : stride何個飛ばしで読み込むか(wavepacket内部のステップを割り切れるかどうかで判定している) (default:1)
-load-min LOAD-MIN : 読み込むファイルの最小値を決めている(wavepacket内部のステップのどこからどこを読み込むか)(default:0)
-load-max LOAD-MAX : 読み込むファイルの最大値を決めている(wavepacket内部のステップのどこからどこを読み込むか)(default:10000000)
-cutoff-au CUTOFF : elses-generate-cubefileで使用するcutoffの値(default:8.0)
-input_mesh_grid PATH : 別階層にあるinput_mesh_gridを現在のディレクトリに移動する(すでにある場合は何もしない)(自動でinput_mesh_gridを作製するオプションと併用した場合はこのオプションは動作しない)
-alpha ALPHA : 最初のステップのmeanと最後のステップの合計のMSDを使用してmean +- alpha*sqrt(MSD)の範囲(正方形)のinput_mesh_gridを作製(既にあるinput_mesh_gridを上書してしまうので注意)
(--position と -alphaを同時に使用すると最大でpositionで作製される範囲になるようにalphaで作製される)
-mesh MESH : オプションでinput_mesh_gridを作製する際のメッシュ数を変更する(メッシュ数は大体 長さ[a.u.]*MESH で作製される) default : 1.0
--parallel : 読み込みとcharのcubeファイル作製をpythonで並列化させる(メモリが十分にあるときだけ使用すること)(並列数は自動でOMP_NUM_THREADSの値になります。OMP_NUM_THREADSがない場合は最大の並列数になります)
-target target : 作りたいwavepacket内部のステップを指定する(-load-min,maxは無視されます)(複数ステップを入力する場合は「115,223,587」のように「,」区切りでスペースがないように入力)
--position : xyzファイルを読み込みその原子が全て含まれるように自動でinput_mesh_gridを作製(既にあるinput_mesh_gridを上書してしまうので注意)
--big-endian : oakleaf,Kでのbynariを読み込む時に使用
--LCAO : LCAO係数データファイル(output_wavefunction)のみを出力するモード(mpi並列で実行しても意味ありません)
-periodic periodic : 周期境界条件を適用するときに使用する。入力方法は -periodic zxy のように拡張したい方向を入力する(順不同)。周期セルに合わせるようにinput_mesh_grid.txtを作成する。ただし、他の自動作成のオプションを使用した場合はそちらが優先される。

出力ファイル名
psi_char_00000000はWaveToolでの通し番号、
step00000000はwavepacket内部のステップ数、
t0000095.79fsは時刻[fs]。


一度作製したoutput_wavefunction.txtを利用してcubeファイルを作成
usage: Main_AutoWave_not_LCAO_20160518.py [-h] [-core CORE] [--parallel]
                                          [-cutoff-au CUTOFF]
                                          <elses-generate-cubefile>
                                          <JSON.CUBE>
例)
python Main_AutoWave_not_LCAO_20160525.py elses-generate-cubefile out.json.cube --parallel
mpi4pyをインストールしていてmpiで並列するとき
mpiexec -n 2 python Main_AutoWave_not_LCAO_20160525.py elses-generate-cubefile out.json.cube --parallel

out.json.cube/real/
out.json.cube/imag/
となっていることが前提。
オプションの意味は、上と同様。
