#ifndef __time_perform_h__
#define __time_perform_h__


#define MAX_TIME_IDX		    128
#define TIME_START_IDX			0
#define TIME_FINAL_IDX			37

#define TOTAL_EXEC_TIME				0
#define READ_AND_INIT_TIME			1
#define WRIT_AND_FINA_TIME			2
#define BOUND_POINT_INST_TIME		3
#define INNER_POINT_INST_TIME	    4
#define LSSP_RECV_EDGE_TIME			5
#define NOSP_RECV_EDGE_TIME			6
#define NOSP_RECV_FACE_TIME			7
#define CONF_RECV_EDGE_TIME			8
#define CONF_RECV_FACE_TIME			9
#define CONS_RECV_EDGE_TIME			10
#define CONS_RECV_FACE_TIME			11
#define CONS_RECV_PNT_RMV_TIME		12
#define QUAL_IMPRV_TIME				13
#define TYPE_DEL_ELEM_TIME			14
#define QUAL_RECV_SH_TRANS_TIME		15
#define MESH_SMOOTHING_TIME			16
#define EDGE_CONTRACT_TIME			17
#define EDGE_SPLITTING_TIME			18
#define QUALITY_EVAL_TIME			19
#define FORM_INIT_CAVITY_TIME		20
#define MODIFY_CAVITY_TIME			21
#define UPDATE_TRIANGULATION_TIME	22
#define CREATE_INNER_POINT_TIME		23
#define CREATE_ELEM_PARA_TIME		24
#define ADD_SUCC_CLEAN_TIME			25
#define ADD_FAIL_CLEAN_TIME			26
#define SPR_CALLING_TIME			27

#define RECV_EDG_LOOP_TIME			28
#define RECV_FAC_LOOP_TIME			29
#define REMESH_ASHELL_TIME			30
#define REMOV_AN_EDGE_TIME			31
#define REFLEX_EDG_CH_TIME			32
#define FIND_ALLSHELL_TIME			33
#define DYN_RSHL_MAIN_TIME			34
#define CHK_EDGE_RECV_TIME			35
#define CHK_FACE_RECV_TIME			36
#define SHELL_TRANS_TIME			37


extern double time_start[MAX_TIME_IDX];
extern double time_final[MAX_TIME_IDX];
extern double time_elaps[MAX_TIME_IDX];
extern int sphere_size_eval;

#endif /* __time_perform_h__ */