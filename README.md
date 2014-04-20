MolecularSimulationTool
====

Currently, followings are implemented.
*  MoleculearDynamics (MD) 
  * QM-MD(Molpro only), MM-MD, QM/MM-MD
  * Excited state Simulation is also OK. Tully's Surface Hopping Algorithm is adopted. 
*  MonteCarlo (MC)
  * Metropolis method is adopted. 

In the future, following functions will be implemented.
* The Format for the other famous quantum chemistry calculation packaged like GAUSSIAN or GAMESS 
* Visualization Tool      

研究用の自作分子シミュレーションプログラムです。
原子間に互いに作用する力を計算しながら、原子群の運動の時間発展を追う分子動力学(MD)法と、乱数を発生させ、熱力学的平衡状態を再現するモンテカルロ(MC)法を実装しています。
全体はpythonで記述していますが、律速箇所である相互作用エネルギーの計算をFortranで書くことで実行速度の改善を行っています。
現段階は、adhocなプログラムの集積なので、今後はより汎用的な用途に対応していくことを目的としています。

srcディクレクトリ内で python run.pyを行うとプログラムが走ります。例ではar結晶のモンテカルロ計算となっています。
mc_energy.datに系のエネルギーを、mc.xyzに系の座標を随時出力していきます。

