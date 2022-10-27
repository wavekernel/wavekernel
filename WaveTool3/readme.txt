���s�R�}���h
python Main_AutoWave.py <elses-generate-cubefile> <JSON> <XYZ> <BASIS>
��)
python Main_AutoWave.py elses-generate-cubefile out.json position.xyz output_basis_information.txt -alpha 2.0 -s 10 --parallel --position
mpi4py���C���X�g�[�����Ă���mpi�ŕ��񂷂�Ƃ�
mpiexec -n 2 python Main_AutoWave.py elses-generate-cubefile out.json position.xyz output_basis_information.txt -alpha 2.0 -s 10 --parallel --position

�I�v�V����
-h : help
-s STRIDE : stride����΂��œǂݍ��ނ�(wavepacket�����̃X�e�b�v������؂�邩�ǂ����Ŕ��肵�Ă���) (default:1)
-load-min LOAD-MIN : �ǂݍ��ރt�@�C���̍ŏ��l�����߂Ă���(wavepacket�����̃X�e�b�v�̂ǂ�����ǂ���ǂݍ��ނ�)(default:0)
-load-max LOAD-MAX : �ǂݍ��ރt�@�C���̍ő�l�����߂Ă���(wavepacket�����̃X�e�b�v�̂ǂ�����ǂ���ǂݍ��ނ�)(default:10000000)
-cutoff-au CUTOFF : elses-generate-cubefile�Ŏg�p����cutoff�̒l(default:8.0)
-input_mesh_grid PATH : �ʊK�w�ɂ���input_mesh_grid�����݂̃f�B���N�g���Ɉړ�����(���łɂ���ꍇ�͉������Ȃ�)(������input_mesh_grid���쐻����I�v�V�����ƕ��p�����ꍇ�͂��̃I�v�V�����͓��삵�Ȃ�)
-alpha ALPHA : �ŏ��̃X�e�b�v��mean�ƍŌ�̃X�e�b�v�̍��v��MSD���g�p����mean +- alpha*sqrt(MSD)�͈̔�(�����`)��input_mesh_grid���쐻(���ɂ���input_mesh_grid���㏑���Ă��܂��̂Œ���)
(--position �� -alpha�𓯎��Ɏg�p����ƍő��position�ō쐻�����͈͂ɂȂ�悤��alpha�ō쐻�����)
-mesh MESH : �I�v�V������input_mesh_grid���쐻����ۂ̃��b�V������ύX����(���b�V�����͑�� ����[a.u.]*MESH �ō쐻�����) default : 1.0
--parallel : �ǂݍ��݂�char��cube�t�@�C���쐻��python�ŕ��񉻂�����(���������\���ɂ���Ƃ������g�p���邱��)(���񐔂͎�����OMP_NUM_THREADS�̒l�ɂȂ�܂��BOMP_NUM_THREADS���Ȃ��ꍇ�͍ő�̕��񐔂ɂȂ�܂�)
-target target : ��肽��wavepacket�����̃X�e�b�v���w�肷��(-load-min,max�͖�������܂�)(�����X�e�b�v����͂���ꍇ�́u115,223,587�v�̂悤�Ɂu,�v��؂�ŃX�y�[�X���Ȃ��悤�ɓ���)
--position : xyz�t�@�C����ǂݍ��݂��̌��q���S�Ċ܂܂��悤�Ɏ�����input_mesh_grid���쐻(���ɂ���input_mesh_grid���㏑���Ă��܂��̂Œ���)
--big-endian : oakleaf,K�ł�bynari��ǂݍ��ގ��Ɏg�p
--LCAO : LCAO�W���f�[�^�t�@�C��(output_wavefunction)�݂̂��o�͂��郂�[�h(mpi����Ŏ��s���Ă��Ӗ�����܂���)
-periodic periodic : �������E������K�p����Ƃ��Ɏg�p����B���͕��@�� -periodic zxy �̂悤�Ɋg����������������͂���(���s��)�B�����Z���ɍ��킹��悤��input_mesh_grid.txt���쐬����B�������A���̎����쐬�̃I�v�V�������g�p�����ꍇ�͂����炪�D�悳���B

�o�̓t�@�C����
psi_char_00000000��WaveTool�ł̒ʂ��ԍ��A
step00000000��wavepacket�����̃X�e�b�v���A
t0000095.79fs�͎���[fs]�B


��x�쐻����output_wavefunction.txt�𗘗p����cube�t�@�C�����쐬
usage: Main_AutoWave_not_LCAO_20160518.py [-h] [-core CORE] [--parallel]
                                          [-cutoff-au CUTOFF]
                                          <elses-generate-cubefile>
                                          <JSON.CUBE>
��)
python Main_AutoWave_not_LCAO_20160525.py elses-generate-cubefile out.json.cube --parallel
mpi4py���C���X�g�[�����Ă���mpi�ŕ��񂷂�Ƃ�
mpiexec -n 2 python Main_AutoWave_not_LCAO_20160525.py elses-generate-cubefile out.json.cube --parallel

out.json.cube/real/
out.json.cube/imag/
�ƂȂ��Ă��邱�Ƃ��O��B
�I�v�V�����̈Ӗ��́A��Ɠ��l�B
