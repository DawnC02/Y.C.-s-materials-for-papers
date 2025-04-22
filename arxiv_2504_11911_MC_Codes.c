///This file is the code of MC simulations for https://doi.org/10.48550/arXiv.2504.11911
///Write by Chen Yi-Duo

///Most notes are written in Chinese, you can read it with translator. The file encoding format is UTF-8.
///If I have spare time in the fulture, I will add more notes in English and upload to Github again.

///Write for using it on Linux system. If you want to use it on Windows system, you have to change the headers.
///I write it on Windows PC. I used CRLF (\r\n).

///About output files. Time series in the folder named "time_series". "result_pt.txt" is useless for results in the paper.
///"result_1.txt" has all results. "result_a" has results averaged over 300 repeats (not used in paper).


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>



//// xorshift128 algorithm by Marsaglia
#include <stdint.h>
static uint32_t w = 88675123;
// return an integer in [0, 2^32) = [0, 4294967296)
uint32_t xor128(void) {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    uint32_t t = x ^ (x << 11);
    x = y;    y = z;    z = w;
    return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}
// return an integer in [0, range)
uint32_t xor128_range(uint32_t range){
    uint32_t threshold = -range % range, r;
    for ( ; ; ) {
        r = xor128();
        if (r >= threshold)return r % range;
    }
}
// return a double number in [0, 1)
#define xor128_double()  (xor128() / 4294967296.0)
// seeding the random number generator
void seed_xor128(uint32_t seed) {
    w ^= seed;
    for(uint16_t i=0x7ff; --i; xor128());  /// warm it up
}


//######################################################################################################################






int** function_URH_maker(int nn,int ll,int gg,int** hyperlink);
int** function_hyperlink(int nn,int ll,int gg,int** hyperlink);
int* function_repeat_test(int ll,int gg,int** hyperlink,int* repeat_result);
int function_i_find(int nn,int gg,double judge);
int** function_matrix_hyperlink_to_player(int nn,int ll,int gg,int** hyperlink);
int** function_matrix_hyperlink_to_player_uni(int nn,int ll,int** hyperlink);
int** function_hyper_graph_maker_2_3(int nn,int ll_2,int ll_3,int** hyperlink,double delta);
int** function_matrix_neighbours_initialization(int nn,int ll,int** hyperlink,int** player);
double function_payoff_calculator(int player_n,int** hyperlink,int** player,double** payoff_2,double** payoff_3,int* strategies,int* game_states,double lower_payoff_delta);
int function_time_step(int nn,int ll,int** hyperlink,int** player,int** neighbours,double** payoff_2,double** payoff_3,int* strategies,double omega,double xct,double transition_probability_repair,double transition_probability_destroy,double lower_payoff_delta,int* game_states);
long int function_count_c(int nn,int* strategies);
int function_folder_create();
int* function_strategies_initialization(int nn,double xc0);
char* function_file_name_txt_density(int xx);
char* function_file_name_txt_time_series(int group_n,int repeat_n);
char* function_file_name_txt_hypergraph_(int group_n,int repeat_n);
int* function_initialize_game_states(int ll);
long int function_count_state_0(int ll,int* game_states);
int function_game_state_update(int player_n,int ll,int nn,int* strategies,int* game_states,int** player,int** hyperlink,double transition_probability_repair,double transition_probability_destroy);
double** function_iteration(int ll,double xc0,int xc0_random_mode,int tt,int nn,double kk,double** payoff_2,double** payoff_3,double delta,double omega,int repeat_n,int group_n,double* zeta_d,int ta,int ll_2,int ll_3,double transition_probability_repair,double transition_probability_destroy,double lower_payoff_delta);
double*** function_main_calculate(double xc0,int xc0_random_mode,int tt,int nn,double kk,double** payoff_2,double** payoff_3,double delta,double omega,int repeat,int group,double group_start,double group_change,int ta,double transition_probability_repair,double transition_probability_destroy,double lower_payoff_delta);


///Do not change anything in this part unless you know its usage! Start.
int zt_pt_mode=0;//分析相变所需统计量功能,0关闭，1开启
int pattern_evo_mode=0;//输出过程中分布情况功能，0关闭，1开启
int transition_c_players=2;//高收益状态需要的合作者数目
///Do not change anything in this part unless you know its usage! End.

int group_calculate_mode=0;//c程序内按组批量计算使用
int divide_sign=47;//Windows: 反斜杠92; others: 正斜杠47


int main() {
    printf("Hello!\n");
    time_t time_start=time(NULL);

    //博弈模拟参数
    int nn=5000;//超图节点数
    int kk=15;//超图节点平均度期望
    double delta=0.6;//3超边占比
    int tt=20000;//总时长
    int ta=5000;//平均量计算时长，即在tt中最后ta时间段计算平均量
    double omega=1./((double)kk);//费米概率选择强度
    int repeat=300;//全同参数重复模拟次数
    int xc0_random_mode=1;//1时随机选取起始合作者密度，0时初始合作者密度为下一个变量
    double xc0=0.5;//初始合作者密度，上一个变量为0时启用
    double lambda=1.;//暂时无用
    double rou=1.;//暂时无用
    double lower_payoff_delta=0.3;//模型收益矩阵参数中的e

    //收益矩阵元素
    double payoff_a=0.8;
    double payoff_b=1.1;
    double payoff_c=0.1;

    if(1){//可以使用"input_p.txt"修改参数，用于批量计算，若txt不存在则不进行修改！！！
        FILE* f_input_parameter=fopen("input_p.txt","r");
        if(f_input_parameter!=NULL){
            fscanf(f_input_parameter,"%lf,%lf,%lf,%lf,%lf",&lower_payoff_delta,&delta,&payoff_a,&payoff_b,&payoff_c);
        }
        fclose(f_input_parameter);
    }

    double payoff_R=payoff_b-payoff_c;
    double payoff_S=-payoff_c;
    double payoff_T=payoff_b;
    double payoff_P=0.;
    double payoff_G=payoff_a/4.+0.5;
    double payoff_W=-payoff_a/4.+0.5;

    //2超边与3超边占比计算
    int ll_2,ll_3;
    if(delta==1.0){
        ll_2=0;
        ll_3=(int)(nn*kk/3.);
    }else if(delta==0.0){
        ll_3=0;
        ll_2=(int)(nn*kk/2.);
    }else{
        ll_3=(int)(nn*kk/(2.*(1.-delta)/delta+3));
        ll_2=(int)(ll_3*(1.-delta)/delta);
    }
    //printf("%d,%d\n",ll_2,ll_3);
    int ll=ll_3+ll_2;

    if(ta>tt){
        printf("Wrong!\nta>tt!!!\n");
        getchar();
        getchar();
        return 0;
    }

    int group=1;//可以用于连续改变单参数，不需要时取“1”！！！
    double group_start=0.05;
    double group_change=0.05;



    double transition_probability_repair=lambda;//修复超边环境的概率
    double transition_probability_destroy=rou;//破坏超边环境的概率
    ///not used

    //收益矩阵
    double** payoff_2=(double**) malloc(2*sizeof(double*));
    double** payoff_3=(double**) malloc(2*sizeof(double*));
    payoff_2[0]=(double*)malloc(2*sizeof(double));
    payoff_2[1]=(double*)malloc(2*sizeof(double));
    payoff_3[0]=(double*)malloc(3*sizeof(double));
    payoff_3[1]=(double*)malloc(3*sizeof(double));
    //双人博弈收益矩阵，程序内定义有所修改，π ij指标i为focal player策略，指标j为对手中合作者数目
    payoff_2[0][0]=payoff_S;//CD
    payoff_2[0][1]=payoff_R;//CC
    payoff_2[1][0]=payoff_P;//DD
    payoff_2[1][1]=payoff_T;//DC
    //三人博弈收益矩阵
    payoff_3[0][0]=payoff_S;//C vs 2D
    payoff_3[0][1]=payoff_G;//C vs C+D
    payoff_3[0][2]=payoff_R;//C vs 2C
    payoff_3[1][0]=payoff_P;//D vs 2D
    payoff_3[1][1]=payoff_W;//D vs C+D
    payoff_3[1][2]=payoff_T;//D vs 2C

    FILE* fpa=fopen("parameter.txt","w");
    fprintf(fpa,"ll=%d,xc0=%f,xc0_random=%d,tt=%d,kk=%d,delta=%f,omega=%f,ta=%d,nn=%d,ll_2=%d,ll_3=%d,lambda=%f,rou=%f,ee=%f\n",ll,xc0,xc0_random_mode,tt,kk,delta,omega,ta,nn,ll_2,ll_3,lambda,rou,lower_payoff_delta);
    fprintf(fpa,"Payoff:R,S,T,P,G,W,a=%f,%f,%f,%f,%f,%f,%f\n",payoff_R,payoff_S,payoff_T,payoff_P,payoff_G,payoff_W,payoff_a);
    fprintf(fpa,"repeat=%d,group=%d,group_start=%f,group_change=%f\n",repeat,group,group_start,group_change);
    fprintf(fpa,"group_cal_mode=%d,state_update_prob_mode=%d,nonlinear_interaction=%d\n",group_calculate_mode,update_state_with_probability_mode,nonlinear_interaction_destroy);
    fprintf(fpa,"使用c程序内group变换delta时此处记录的ll与l2与l3无效！\n");

    function_folder_create();//创建存储数据需要的文件夹

    double*** final_result= function_main_calculate(xc0,xc0_random_mode,tt,nn,kk,payoff_2,payoff_3,delta,omega,repeat,group,group_start,group_change,ta,transition_probability_repair,transition_probability_destroy,lower_payoff_delta);

    FILE* fxn= fopen("result_1.txt","w");
    for(int g=0;g<group;g++){
        for(int r=0;r<repeat;r++){
            fprintf(fxn,"%d\t%d\t%f\t%f\t%f\t%f\n",g+1,r+1,final_result[r][g][0],final_result[r][g][1],final_result[r][g][2],final_result[r][g][3]);
            //fprintf(fxn,"%d\t%d\t%f\t%f\t%f\n",g+1,r+1,final_result[r][g][0],final_result[r][g][1],final_result[r][g][2]);
        }
    }
    fclose(fxn);

    FILE* fxa= fopen("result_a.txt","w");
    for(int g=0;g<group;g++){
        double xcaa=0.,naa=0.;
        for(int r=0;r<repeat;r++){
            xcaa+=final_result[r][g][2];
            naa+=final_result[r][g][3];
        }
        if(group_calculate_mode){
            fprintf(fxa,"%d\t%f\t%f\t%f\n",g+1,group_start+g*group_change,xcaa/repeat,naa/repeat);
        }else{
            fprintf(fxa,"%d\t%f\t%f\t%f\n",g+1,delta,xcaa/repeat,naa/repeat);
        }

    }

    for(double*** pi=final_result;pi<final_result+repeat;pi++){
        for(double** pj=*pi;pj<*pi+group;pj++){
            free(*pj);
            *pj=NULL;
        }
        free(*pi);
        *pi=NULL;
    }
    free(final_result);
    final_result=NULL;

    printf("Finished!\nTime spent %ds!\nTime: %d\n",(int)(time(NULL)-time_start),(int)time(NULL));
    fprintf(fpa,"\nTime spent %ds!\n",(int)(time(NULL)-time_start));
    fclose(fpa);

    free(payoff_2[0]);
    free(payoff_2[1]);
    free(payoff_3[0]);
    free(payoff_3[1]);
    free(payoff_2);
    free(payoff_3);

    return 0;
}


double*** function_main_calculate(double xc0,int xc0_random_mode,int tt,int nn,double kk,double** payoff_2,double** payoff_3,double delta,double omega,int repeat,int group,double group_start,double group_change,int ta,double transition_probability_repair,double transition_probability_destroy,double lower_payoff_delta){
    double*** final_result=(double***)malloc(repeat*sizeof(double**));
    //最终结果文件，访问方式 r g （gs xc0 xc n）
    for(double*** pi=final_result;pi<final_result+repeat;pi++){
        *pi=(double**)malloc(group*sizeof(double*));
        for(double** pj=*pi;pj<*pi+group;pj++){
            *pj=(double*)malloc(4*sizeof(double));
        }
    }

    double** result_pt=(double**)malloc(group*sizeof(double*));
    for(double** p_rpt_i=result_pt;p_rpt_i<result_pt+group;p_rpt_i++){
        *p_rpt_i=(double*)malloc(5*sizeof(double));
    }

    //FILE* fb= fopen("result_bif.txt","w");//可以绘制长时间多系统分岔图，需要将所有相关代码取消注释

    for(int g=0;g<group;g++){
        double pt_m0=0.,pt_chi0=0.,pt_u0=0.;

        if(group_calculate_mode){
            delta=group_start+g*group_change;//使用group变化变量，启用时取消注释本行代码
        }
        for(int r=0;r<repeat;r++){
            //2超边与3超边占比计算
            int ll_2,ll_3;
            if(delta==1.0){
                ll_2=0;
                ll_3=(int)(nn*kk/3.);
            }else if(delta==0.0){
                ll_3=0;
                ll_2=(int)(nn*kk/2.);
            }else{
                ll_3=(int)(nn*kk/(2.*(1.-delta)/delta+3));
                ll_2=(int)(ll_3*(1.-delta)/delta);
            }
            //printf("%d,%d\n",ll_2,ll_3);
            int ll=ll_3+ll_2;

            double* zeta_d=(double*)malloc(3*sizeof(double));
            for(int zt_i=0;zt_i<3;zt_i++){
                zeta_d[zt_i]=0.;
            }

            if(xc0_random_mode){
                xc0=xor128_double();
            }
            //printf("%f\n",delta);
            double** time_series= function_iteration(ll,xc0,xc0_random_mode,tt,nn,kk,payoff_2,payoff_3,delta,omega,r,g,zeta_d,ta,ll_2,ll_3,transition_probability_repair,transition_probability_destroy,lower_payoff_delta);

            if(zt_pt_mode){
                pt_m0+=zeta_d[0];
                pt_chi0+=zeta_d[1]-zeta_d[0]*zeta_d[0];
                pt_u0+=zeta_d[2]/(3.*zeta_d[1]*zeta_d[1]);
            }

            free(zeta_d);
            zeta_d=NULL;


            double xa=0.,na=0.;
            char* tsname=function_file_name_txt_time_series(g,r);
            FILE* fp= fopen(tsname,"w");
            free(tsname);
            tsname=NULL;
            for(int t=0;t<tt+1;t++){
                fprintf(fp,"%d\t%f\t%f\n",t,time_series[t][0],time_series[t][1]);
                //fprintf(fp,"%d\t%f\n",t,time_series[t][0]);
                if(t>tt-ta){
                    xa+=time_series[t][0]/ta;
                    na+=time_series[t][1]/ta;
                }
                /*if(t>tt-500){
                    fprintf(fb,"%f\t%f\t%f\n",mu,time_series[t][0],time_series[t][1]);
                }*/
            }
            fclose(fp);
            for(double** ppi=time_series;ppi<time_series+tt+1;ppi++){
                free(*ppi);
                *ppi=NULL;
            }
            free(time_series);
            time_series=NULL;

            if(group_calculate_mode){
                final_result[r][g][0]=group_start+g*group_change;
            }else{
                final_result[r][g][0]=delta;//
            }
            final_result[r][g][1]=xc0;
            final_result[r][g][2]=xa;
            final_result[r][g][3]=na;
        }

        if(zt_pt_mode){
            pt_m0/=(double)repeat;
            pt_chi0/=(double)repeat;
            pt_u0/=(double)repeat;

            result_pt[g][0]=(double)nn;
            result_pt[g][1]=delta;
            result_pt[g][2]=pt_m0;
            result_pt[g][3]=nn*pt_chi0;
            result_pt[g][4]=1.-pt_u0;
        }
    }

    if(zt_pt_mode){
        FILE* f_rpt= fopen("result_pt.txt","w");
        for(int fp_g=0;fp_g<group;fp_g++){
            fprintf(f_rpt,"%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",result_pt[fp_g][0],result_pt[fp_g][1],result_pt[fp_g][2],result_pt[fp_g][3],result_pt[fp_g][4]);
        }
        fclose(f_rpt);
    }

    //fclose(fb);
    for(double** pfpt=result_pt;pfpt<result_pt+group;pfpt++){
        free(*pfpt);
        *pfpt=NULL;
    }
    free(result_pt);
    result_pt=NULL;

    return final_result;
}


double** function_iteration(int ll,double xc0,int xc0_random_mode,int tt,int nn,double kk,double** payoff_2,double** payoff_3,double delta,double omega,int repeat_n,int group_n,double* zeta_d,int ta,int ll_2,int ll_3,double transition_probability_repair,double transition_probability_destroy,double lower_payoff_delta){
    int pro_time=(tt>1000)?(tt/100):(tt/10);

    time_t time_start_group=time(NULL);

    int** hyperlink=(int**)malloc((ll_2+ll_3)*sizeof(int*));
    //hyperlink矩阵，存储超图内超边信息，访问指标ij，i指标为超边编号，j=0为超边内个体数目，j=1，2，……为超边内个体编号
    int middle_h=0;
    for(int** pi=hyperlink;pi<hyperlink+ll_2+ll_3;pi++){
        if(middle_h<ll_2){
            *pi=(int*)malloc(3*sizeof(int));
        }else{
            *pi=(int*)malloc(4*sizeof(int));
        }
        middle_h++;
    }

    int** player= function_hyper_graph_maker_2_3(nn,ll_2,ll_3,hyperlink,delta);
    //player矩阵，存储超图内个体信息，访问指标ij，i指标为个体编号，j=0为个体节点的度，j=1，2，……为个体所在超边编号

    int** neighbours= function_matrix_neighbours_initialization(nn,ll,hyperlink,player);
    //neighbours矩阵，存储超图内个体邻居信息，访问指标ij，i指标为个体编号，j=0为个体邻居总数，j=1，2，……为个体邻居的编号

    int* strategies;
    //个体策略态矢量，存储个体的策略信息，访问指标为个体编号
    strategies= function_strategies_initialization(nn,xc0);

    int* game_states= function_initialize_game_states(ll);

    double** time_series=(double**)malloc((tt+1)*sizeof(double*));
    for(double** pi=time_series;pi<time_series+tt+1;pi++){
        *pi=(double*)malloc(2*sizeof(double));//可以增加其他统计量，预留
    }
    double** pt=time_series;

    **pt= (double)function_count_c(nn,strategies)/nn;
    *(*pt+1)=(double)function_count_state_0(ll,game_states)/ll;
    pt++;


    if(0){
        FILE* fp00= fopen("distribution_start.txt","w");
        int* pi=strategies;
        //double** pei=system_env;
        for(int x=0;x<nn;x++){

            fprintf(fp00,"%d\t%d\n",x,*pi);

            pi++;
            //pei++;
        }
        fclose(fp00);
    }//输出开始时的分布情况，多组重复时会相互覆盖

    double xct0=**time_series;

    for(int t=0;t<tt;t++){
        if(t%500==0){
            seed_xor128(time(NULL)+ xor128_range(1000));
        }

        for(int i=0;i<nn;i++){
            function_time_step(nn,ll,hyperlink,player,neighbours,payoff_2,payoff_3,strategies,omega,xct0,transition_probability_repair,transition_probability_destroy,lower_payoff_delta,game_states);
        }

        double xc= (double)function_count_c(nn,strategies)/nn;
        double env= (double)function_count_state_0(ll,game_states)/ll;
        xct0=xc;
        //printf("%f\n",xc);

        **(pt)=xc;
        *(*pt+1)=env;
        pt++;

        if((t+1)%pro_time==0){//进度输出
            FILE* fmp= fopen("progress.txt","w");
            fprintf(fmp,"group %d:\trepeat %d:\t%.3f percent done!\tTime: %ds!\n",group_n+1,repeat_n+1,(double)(t+1)*100./tt,(int)(time(NULL)-time_start_group));
            fclose(fmp);
        }

        if((zt_pt_mode)&&(t>=tt-ta)){
            double zeta=xc;
            zeta_d[0]+=zeta;
            zeta_d[1]+=zeta*zeta;
            zeta_d[2]+=zeta*zeta*zeta*zeta;
        }

        if(pattern_evo_mode){
            //不使用时本代码块从此开始注掉，仅无重复时有效使用
            if((((t+1)%5==0)&&(t+1<=tt)) || (((t+1)%1==0)&&(t+1>=tt+200))){
                //记录不同时间下系统内物种分布情况，仅在无重复时使用
                char* file_name_m=function_file_name_txt_density(t+1);
                FILE* fp0= fopen(file_name_m,"w");
                free(file_name_m);

                int* pi=strategies;
                //double** pei=system_env;
                for(int x=0;x<ll;x++){

                    //fprintf(fp0,"%d\t%d\t%d\t%f\n",x,y,*pi,*pei);
                    fprintf(fp0,"%d\t%d\n",x,*pi);

                    pi++;
                    //pei++;
                }

                fclose(fp0);

            }//代码块结束*/
        }

    }

    if(0){
        //输出该组的超图结构

        char* fcname= function_file_name_txt_hypergraph_(group_n,repeat_n);
        FILE* fpo_h= fopen(fcname,"w");
        free(fcname);
        fcname=NULL;

        int pom=0;
        for(int** ppp=hyperlink;ppp<hyperlink+ll_3+ll_2;ppp++){
            fprintf(fpo_h,"%d\t%d",pom,**ppp);
            for(int* ppp1=*ppp+1;ppp1<*ppp+1+**ppp;ppp1++){
                fprintf(fpo_h,"\t%d",*ppp1);
            }
            fprintf(fpo_h,"\n");
            pom++;
        }

        fclose(fpo_h);
    }

    for(int** pi=hyperlink;pi<hyperlink+ll_2+ll_3;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(hyperlink);
    hyperlink=NULL;

    for(int** pi=player;pi<player+nn;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(player);
    player=NULL;

    for(int** pi=neighbours;pi<neighbours+nn;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(neighbours);
    neighbours=NULL;

    free(strategies);
    strategies=NULL;

    for(int i=0;i<3;i++){
        zeta_d[i]/=ta;
    }

    free(game_states);
    game_states=NULL;

    return time_series;
}


int function_time_step(int nn,int ll,int** hyperlink,int** player,int** neighbours,double** payoff_2,double** payoff_3,int* strategies,double omega,double xct,double transition_probability_repair,double transition_probability_destroy,double lower_payoff_delta,int* game_states){
    int focal_player=(int) xor128_range(nn);
    //if((xct>(double)(nn-1.)/nn)||(xct<(double)1./nn)){}

    if(xor128_double()<(double)(2./nn)){
        //  2/N噪声，随机选择策略
        strategies[focal_player]= (int)xor128_range(2);
        return 0;
    }

    int model_player=neighbours[focal_player][(int) xor128_range(neighbours[focal_player][0])+1];
    //model player为focal player邻居中随机选择

    double payoff_f= function_payoff_calculator(focal_player,hyperlink,player,payoff_2,payoff_3,strategies,game_states,lower_payoff_delta);
    double payoff_m= function_payoff_calculator(model_player,hyperlink,player,payoff_2,payoff_3,strategies,game_states,lower_payoff_delta);

    function_game_state_update(focal_player,ll,nn,strategies,game_states,player,hyperlink,transition_probability_repair,transition_probability_destroy);
    function_game_state_update(model_player,ll,nn,strategies,game_states,player,hyperlink,transition_probability_repair,transition_probability_destroy);

    double update_probability=1./(1.+exp(omega*(payoff_f-payoff_m)));
    double update_judge=xor128_double();

    //printf("%d,%d,%d,%d,%f,%f,%f,%f\n",focal_player,model_player,strategies[focal_player],strategies[model_player],payoff_f,payoff_m,update_probability,update_judge);//debug

    if(update_judge<update_probability){
        strategies[focal_player]=strategies[model_player];
    }

    return 0;
}


int function_game_state_update(int player_n,int ll,int nn,int* strategies,int* game_states,int** player,int** hyperlink,double transition_probability_repair,double transition_probability_destroy){

    for(int* pj=*(player+player_n)+1;pj<*(player+player_n)+1+**(player+player_n);pj++){
        int sop=0;
        for(int* psk=*(hyperlink+*pj)+1;psk<*(hyperlink+*pj)+1+**(hyperlink+*pj);psk++){
            sop+=strategies[*psk];
        }

        int game_state_0=game_states[*pj];

        sop=**(hyperlink+*pj)-sop;//超边内合作者数目
        int k_m=**(hyperlink+*pj);

        if((sop>=transition_c_players) && (game_state_0!=0)){
            //合作者达到一定数量后到达高收益状态
            game_states[*pj]=0;
        }else if((sop<transition_c_players) && (game_state_0==0)){
            //合作者数目不足时到达低收益状态
            game_states[*pj]=1;
        }

    }

    return 0;
}


int* function_initialize_game_states(int ll){
    int* game_states=(int*)malloc(ll*sizeof(int));
    for(int* pi=game_states;pi<game_states+ll;pi++){
        *pi=0;
    }

    return game_states;
}


long int function_count_c(int nn,int* strategies){
    long int count=0;
    for(int* pi=strategies;pi<strategies+nn;pi++){
        if(*pi==0){
            count+=1;
        }
    }
    return count;
}


long int function_count_state_0(int ll,int* game_states){
    long int count=0;
    for(int* pi=game_states;pi<game_states+ll;pi++){
        if(*pi==0){
            count+=1;
        }
    }
    return count;
}


int* function_strategies_initialization(int nn,double xc0){
    int* strategies=(int*)malloc(nn*sizeof(int));
    for(int* pi=strategies;pi<strategies+nn;pi++){
        *pi=(xor128_double()<xc0)?(0):(1);
    }

    return strategies;
}


double function_payoff_calculator(int player_n,int** hyperlink,int** player,double** payoff_2,double** payoff_3,int* strategies,int* game_states,double lower_payoff_delta){
    int s0=strategies[player_n];
    int c0=1-s0;
    double payoff0=0.;
    for(int* pj=*(player+player_n)+1;pj<*(player+player_n)+1+**(player+player_n);pj++){
        int sop=0;
        for(int* psk=*(hyperlink+*pj)+1;psk<*(hyperlink+*pj)+1+**(hyperlink+*pj);psk++){
            sop+=strategies[*psk];
        }

        int game_state_0=game_states[*pj];

        sop=**(hyperlink+*pj)-sop;//超边内合作者数目
        sop-=c0;//减去自己

        if(**(hyperlink+*pj)==2){
            payoff0+=payoff_2[s0][sop];
        }else{
            payoff0+=payoff_3[s0][sop];
        }

        if((sop!=0)&&game_state_0==1){
            payoff0-=lower_payoff_delta;
        }
    }

    return payoff0;
}


int** function_matrix_neighbours_initialization(int nn,int ll,int** hyperlink,int** player){
    //统计生成每个玩家的邻居数目与具体邻居矩阵
    //因为有3阶或以上超边存在，可能存在单个邻居被反复连接，按超边选取选择概率不均匀

    int** neighbours=(int**)malloc(nn*sizeof(int*));
    int** ppi=player;
    for(int** pi=neighbours;pi<neighbours+nn;pi++){
        *pi=(int*)malloc((2*(**ppi)+1)*sizeof(int));
        **pi=0;
        ppi++;
    }

    int i=0;
    ppi=player;
    for(int** pi=neighbours;pi<neighbours+nn;pi++){
        int* pj=*pi+1;
        int neighbour_count=0;
        for(int* ppj=*ppi+1;ppj<*ppi+1+**ppi;ppj++){
            for(int* pek=*(hyperlink+*ppj)+1;pek<*(hyperlink+*ppj)+1+**(hyperlink+*ppj);pek++){
                if(*pek!=i){
                    int judge_not_in=1;
                    for(int* ptk=*pi+1;ptk<*pi+neighbour_count+1;ptk++){
                        if(*ptk==*pek){
                            judge_not_in=0;
                            break;
                        }
                    }

                    if(judge_not_in){
                        *pj++=*pek;
                        neighbour_count++;
                    }
                }

            }
        }

        ppi++;
        i++;
        **pi=neighbour_count;
    }

    return neighbours;
}


int** function_hyper_graph_maker_2_3(int nn,int ll_2,int ll_3,int** hyperlink,double delta){
    //生成2/3超阶混合随机超图
    //注意：注掉的代码块需要放在该函数调用前！！！
    /*

    int** hyperlink=(int**)malloc((ll_2+ll_3)*sizeof(int*));
    int middle_h=0;
    for(int** pi=hyperlink;pi<hyperlink+ll_2+ll_3;pi++){
        if(middle_h<ll_2){
            *pi=(int*)malloc(3*sizeof(int));
        }else{
            *pi=(int*)malloc(4*sizeof(int));
        }
        middle_h++;
    }

     */

    int** hyperlink_2=(int**)malloc(ll_2* sizeof(int*));//存储每个连接涉及的个体
    for(int** pi=hyperlink_2;pi<hyperlink_2+ll_2;pi++){
        *pi=(int*)malloc(2*sizeof(int));
    }
    int** hyperlink_3=(int**)malloc(ll_3* sizeof(int*));//存储每个连接涉及的个体
    for(int** pi=hyperlink_3;pi<hyperlink_3+ll_3;pi++){
        *pi=(int*)malloc(3*sizeof(int));
    }


    //使用近似算法生成均匀随机超图后，需要考虑孤立节点与孤立团簇问题，均使用单边打断重连解决。

    int** player_2= function_URH_maker(nn,ll_2,2,hyperlink_2);
    int** player_3= function_URH_maker(nn,ll_3,3,hyperlink_3);

    //首先考虑孤立节点问题

    for(int i=0;i<nn;i++){

        if(player_2[i][0]==0){
            if(player_3[i][0]==0){
                //孤立节点

                //printf("Not connected node!\n");

                int judge_replaced=1;
                int replaced;
                while(judge_replaced){
                    replaced= (int)xor128_range(nn);
                    if(player_2[replaced][0]+player_3[replaced][0]>=2){
                        judge_replaced=0;
                    }
                }

                int middle_k2=player_2[replaced][0],middle_k3=player_3[replaced][0];
                int judge_change_replaced_edge= (int)xor128_range(middle_k2+middle_k3);
                if(judge_change_replaced_edge<middle_k2) {
                    int replaced_hyperedge = player_2[replaced][judge_change_replaced_edge+1];
                    for(int* ri=*(hyperlink_2+replaced_hyperedge);ri<*(hyperlink_2+replaced_hyperedge)+2;ri++){
                        if(*ri==replaced){
                            *ri=i;
                        }
                    }

                    for(int** fpi=player_2;fpi<player_2+nn;fpi++){
                        free(*fpi);
                        *fpi=NULL;
                    }
                    free(player_2);
                    player_2=NULL;

                    player_2= function_matrix_hyperlink_to_player(nn,ll_2,2,hyperlink_2);
                }else{
                    int replaced_hyperedge = player_3[replaced][judge_change_replaced_edge-middle_k2+1];
                    for(int* ri=*(hyperlink_3+replaced_hyperedge);ri<*(hyperlink_3+replaced_hyperedge)+3;ri++){
                        if(*ri==replaced){
                            *ri=i;
                        }
                    }

                    for(int** fpi=player_3;fpi<player_3+nn;fpi++){
                        free(*fpi);
                        *fpi=NULL;
                    }
                    free(player_3);
                    player_3=NULL;

                    player_3= function_matrix_hyperlink_to_player(nn,ll_3,3,hyperlink_3);
                }

            }
        }
    }

    if(*player_2==NULL){
        printf("*player_2==NULL\n");//debug
    }
    if(*player_3==NULL){
        printf("*player_3==NULL\n");//debug
    }


    //然后考虑孤立团簇问题，发现非全覆盖团簇时，随机选择一个团簇外节点和一个团簇内节点建立一个两个体超边，打破一个不影响连通性的两个体超边，若ll2=0，使用三个体超边完成上述过程
    int cluster_size=0;
    while(cluster_size!=nn){
        cluster_size=0;
        int* cluster_connected_nodes=(int*)malloc(nn*sizeof(int));
        int finished_nodes=0;
        int* pfn=cluster_connected_nodes;

        cluster_size++;
        *pfn++=0;

        while(finished_nodes<cluster_size){
            int node_middle=cluster_connected_nodes[finished_nodes];
            for(int* pi=*(player_2+node_middle)+1;pi<*(player_2+node_middle)+player_2[node_middle][0]+1;pi++){
                for(int* pj=*(hyperlink_2+*pi);pj<*(hyperlink_2+*pi)+2;pj++){
                    int test_inside=0;

                    for(int* pti=cluster_connected_nodes;pti<cluster_connected_nodes+cluster_size;pti++){
                        if(*pti==*pj){
                            test_inside=1;
                            break;
                        }
                    }

                    if(!test_inside){
                        *pfn++=*pj;
                        cluster_size++;
                    }
                }
            }

            for(int* pi=*(player_3+node_middle)+1;pi<*(player_3+node_middle)+player_3[node_middle][0]+1;pi++){
                for(int* pj=*(hyperlink_3+*pi);pj<*(hyperlink_3+*pi)+3;pj++){
                    int test_inside=0;

                    for(int* pti=cluster_connected_nodes;pti<cluster_connected_nodes+cluster_size;pti++){
                        if(*pti==*pj){
                            test_inside=1;
                            break;
                        }
                    }

                    if(!test_inside){
                        *pfn++=*pj;
                        cluster_size++;
                    }
                }
            }


            finished_nodes++;
        }

        if(cluster_size<nn){
            printf("Cluster size %d smaller than N!\n",cluster_size);
            int judge_not_inside=0;
            int replaced;
            while(!judge_not_inside){
                replaced= (int)xor128_range(nn);
                int test_inside=0;
                for(int* pti=cluster_connected_nodes;pti<cluster_connected_nodes+finished_nodes;pti++){
                    if(*pti==replaced){
                        test_inside=1;
                        break;
                    }
                }

                if(test_inside==0){
                    judge_not_inside=1;
                }
            }

            int inside_connected=cluster_connected_nodes[(int) xor128_range(cluster_size)];

            if(delta<0.9){
                int judge_both_do_not_only_have_one_2_player_edge=0;
                int dying_edge;
                while(!judge_both_do_not_only_have_one_2_player_edge){
                    dying_edge=(int) xor128_range(ll_2);
                    if(player_2[hyperlink_2[dying_edge][0]][0]+player_3[hyperlink_2[dying_edge][0]][0]>=2){
                        if(player_2[hyperlink_2[dying_edge][1]][0]+player_3[hyperlink_2[dying_edge][1]][0]>=2){
                            judge_both_do_not_only_have_one_2_player_edge=1;
                        }
                    }
                }
                hyperlink_2[dying_edge][0]=replaced;
                hyperlink_2[dying_edge][1]=inside_connected;

                for(int** fpi=player_2;fpi<player_2+nn;fpi++){
                    free(*fpi);
                    *fpi=NULL;
                }
                free(player_2);
                player_2=NULL;

                player_2= function_matrix_hyperlink_to_player(nn,ll_2,2,hyperlink_2);
                if(*player_2==NULL){
                    printf("*player_2==NULL\n");//debug
                }
            }else{
                int judge_both_do_not_only_have_one_3_player_edge=0;
                int dying_edge;
                while(!judge_both_do_not_only_have_one_3_player_edge){
                    dying_edge=(int) xor128_range(ll_3);
                    if(player_2[hyperlink_3[dying_edge][0]][0]+player_3[hyperlink_3[dying_edge][0]][0]>=2){
                        if(player_2[hyperlink_3[dying_edge][1]][0]+player_3[hyperlink_3[dying_edge][1]][0]>=2){
                            if(player_2[hyperlink_3[dying_edge][2]][0]+player_3[hyperlink_3[dying_edge][2]][0]>=2){
                                judge_both_do_not_only_have_one_3_player_edge=1;
                            }
                        }
                    }
                }
                int judge_third_player_repeat=1;
                int third_player;
                while(judge_third_player_repeat){
                    third_player=cluster_connected_nodes[(int) xor128_range(cluster_size)];
                    if(third_player!=inside_connected){
                        judge_third_player_repeat=0;
                    }
                }
                hyperlink_3[dying_edge][0]=replaced;
                hyperlink_3[dying_edge][1]=inside_connected;
                hyperlink_3[dying_edge][2]=third_player;

                for(int** fpi=player_3;fpi<player_3+nn;fpi++){
                    free(*fpi);
                    *fpi=NULL;
                }
                free(player_3);
                player_3=NULL;

                player_3= function_matrix_hyperlink_to_player(nn,ll_3,3,hyperlink_3);
                if(*player_3==NULL){
                    printf("*player_3==NULL\n");//debug
                }
            }

        }

        //printf("cluster_size=%d,nn=%d\n",cluster_size,nn);//debug

        free(cluster_connected_nodes);
        cluster_connected_nodes=NULL;

    }


    //重新设计存储超边和玩家信息的矩阵
    int middle_t0=0;
    int** th2=hyperlink_2;
    int** th3=hyperlink_3;
    for(int** th=hyperlink;th<hyperlink+ll_2+ll_3;th++){
        if(middle_t0<ll_2){
            int* thj=*th;
            int* th2j=*th2;
            *thj++=2;
            *thj++=*th2j++;
            *thj=*th2j;
            th2++;
        }else{
            int* thj=*th;
            int* th3j=*th3;
            *thj++=3;
            *thj++=*th3j++;
            *thj++=*th3j++;
            *thj=*th3j;
            th3++;
        }
        middle_t0++;
    }
    int** player= function_matrix_hyperlink_to_player_uni(nn,ll_2+ll_3,hyperlink);

    for(int** pi=hyperlink_2;pi<hyperlink_2+ll_2;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(hyperlink_2);
    hyperlink_2=NULL;

    for(int** pi=player_2;pi<player_2+nn;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(player_2);
    player_2=NULL;

    for(int** pi=hyperlink_3;pi<hyperlink_3+ll_3;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(hyperlink_3);
    hyperlink_3=NULL;

    for(int** pi=player_3;pi<player_3+nn;pi++){
        free(*pi);
        *pi=NULL;
    }
    free(player_3);
    player_3=NULL;

    return player;
}


int** function_matrix_hyperlink_to_player(int nn,int ll,int gg,int** hyperlink){
    //生成URH，首先创建所有超连接让其连接g个个体，然后创建player矩阵并标记每个个体涉及的连接
    //传入hyperlink指针需为NULL
    //player矩阵每个向量为该玩家涉及的超连接数目和超连接代号
    int** player=(int**)malloc(nn*sizeof(int*));
    int* temp_storage=(int*)malloc(nn*sizeof(int));

    for(int* pts=temp_storage;pts<temp_storage+nn;pts++){
        *pts=0;
    }

    for(int** pi=hyperlink;pi<hyperlink+ll;pi++){
        for(int* pj=*pi;pj<*pi+gg;pj++){
            temp_storage[*pj]++;
        }
    }

    /*printf("debug_1\n");
    for(int i=0;i<nn;i++){
        printf("%d\n",temp_storage[i]);
    }*/


    int* ppi0=temp_storage;
    for(int** ppi=player;ppi<player+nn;ppi++,ppi0++){
        *ppi=(int*)malloc((*ppi0+1)*sizeof(int));
    }

    ppi0=temp_storage;
    for(int** ppi=player;ppi<player+nn;ppi++,ppi0++){
        **ppi=*ppi0;
        *ppi0=0;//计数后清空，后面二次利用此数组计数
    }

    int temp_i=0;
    for(int** pi=hyperlink;pi<hyperlink+ll;pi++){
        for(int* pj=*pi;pj<*pi+gg;pj++){
            temp_storage[*pj]++;
            player[*pj][temp_storage[*pj]]=temp_i;
        }
        temp_i++;
    }

    free(temp_storage);
    temp_storage=NULL;

    return player;
}


int** function_matrix_hyperlink_to_player_uni(int nn,int ll,int** hyperlink){
    //生成URH，首先创建所有超连接让其连接g个个体，然后创建player矩阵并标记每个个体涉及的连接
    //传入hyperlink指针需为NULL
    //player矩阵每个向量为该玩家涉及的超连接数目和超连接代号
    int** player=(int**)malloc(nn*sizeof(int*));
    int* temp_storage=(int*)malloc(nn*sizeof(int));

    for(int* pts=temp_storage;pts<temp_storage+nn;pts++){
        *pts=0;
    }

    for(int** pi=hyperlink;pi<hyperlink+ll;pi++){
        int gg=**pi;
        for(int* pj=*(pi)+1;pj<*pi+gg+1;pj++){
            temp_storage[*pj]++;
        }
    }

    /*printf("debug_1\n");
    for(int i=0;i<nn;i++){
        printf("%d\n",temp_storage[i]);
    }*/


    int* ppi0=temp_storage;
    for(int** ppi=player;ppi<player+nn;ppi++,ppi0++){
        *ppi=(int*)malloc((*ppi0+1)*sizeof(int));
    }

    ppi0=temp_storage;
    for(int** ppi=player;ppi<player+nn;ppi++,ppi0++){
        **ppi=*ppi0;
        *ppi0=0;//计数后清空，后面二次利用此数组计数
    }

    int temp_i=0;
    for(int** pi=hyperlink;pi<hyperlink+ll;pi++){
        int gg=**pi;
        for(int* pj=*(pi)+1;pj<*pi+gg+1;pj++){
            temp_storage[*pj]++;
            player[*pj][temp_storage[*pj]]=temp_i;
        }
        temp_i++;
    }

    free(temp_storage);
    temp_storage=NULL;

    return player;
}


int** function_URH_maker(int nn,int ll,int gg,int** hyperlink){
    //生成URH，首先创建所有超连接让其连接g个个体，然后创建player矩阵并标记每个个体涉及的连接
    //传入hyperlink指针需为NULL
    //player矩阵每个向量为该玩家涉及的超连接数目和超连接代号
    int** player=(int**)malloc(nn*sizeof(int*));
    hyperlink= function_hyperlink(nn,ll,gg,hyperlink);
    int* temp_storage=(int*)malloc(nn*sizeof(int));

    for(int* pts=temp_storage;pts<temp_storage+nn;pts++){
        *pts=0;
    }

    for(int** pi=hyperlink;pi<hyperlink+ll;pi++){
        for(int* pj=*pi;pj<*pi+gg;pj++){
            temp_storage[*pj]++;
        }
    }

    /*printf("debug_1\n");
    for(int i=0;i<nn;i++){
        printf("%d\n",temp_storage[i]);
    }*/


    int* ppi0=temp_storage;
    for(int** ppi=player;ppi<player+nn;ppi++,ppi0++){
        *ppi=(int*)malloc((*ppi0+1)*sizeof(int));
    }

    ppi0=temp_storage;
    for(int** ppi=player;ppi<player+nn;ppi++,ppi0++){
        **ppi=*ppi0;
        *ppi0=0;//计数后清空，后面二次利用此数组计数
    }

    int temp_i=0;
    for(int** pi=hyperlink;pi<hyperlink+ll;pi++){
        for(int* pj=*pi;pj<*pi+gg;pj++){
            temp_storage[*pj]++;
            player[*pj][temp_storage[*pj]]=temp_i;
        }
        temp_i++;
    }

    free(temp_storage);
    temp_storage=NULL;

    return player;
}


int** function_hyperlink(int nn,int ll,int gg,int** hyperlink){
    //创建所有超连接让其均连接g个个体
    //生成的超链接里玩家取值为[0,nn)

    seed_xor128(time(NULL)+ xor128_range(1000));

    int** pi=hyperlink;
    for(int i=0;i<ll;i++){

        int* pj=*pi;
        int temp=0;
        for(int j=0;j<gg;j++){

            if(j==0){
                double middle_judge=xor128_double();
                temp+= function_i_find(nn-temp,gg-j,middle_judge);
                *pj=temp;
            }else{
                double middle_judge=xor128_double();
                temp+= function_i_find(nn-temp-1,gg-j,middle_judge)+1;
                *pj=temp;
                if(temp>=nn){
                    printf("wrong_2\n");
                }
            }

            if(j<gg-1){
                pj++;
            }
        }


        if(i<ll-1){
            pi++;
        }
    }

    int repeat_group=1;
    while(repeat_group){

        int* repeat_record= function_repeat_test(ll,gg,hyperlink,&repeat_group);
        if(repeat_group){
            for(int i=0;i<repeat_group;i++){
                int* pj=*(hyperlink+repeat_record[i]);
                int temp=0;
                for(int j=0;j<gg;j++){

                    if(j==0){
                        double middle_judge=xor128_double();
                        temp+= function_i_find(nn-temp,gg-j,middle_judge);
                        *pj=temp;
                    }else{
                        double middle_judge=xor128_double();
                        temp+= function_i_find(nn-temp-1,gg-j,middle_judge)+1;
                        *pj=temp;
                    }

                    if(j<gg-1){
                        pj++;
                    }
                }
            }
        }
        free(repeat_record);
        repeat_record=NULL;
    }

    return hyperlink;
}


int* function_repeat_test(int ll,int gg,int** hyperlink,int* repeat_result){
    //检查重复的超连接
    int* temp_record=(int*)malloc(ll*sizeof(int));
    int repeat_count=0;

    int** p1i=hyperlink;
    int* pr=temp_record;
    for(int i=0;i<ll;i++){
        int** p2i=p1i+1;
        for(int j=i+1;j<ll;j++){
            int repeat_value=0;
            int* p1j=*p1i;
            int* p2j=*p2i;
            for(int k=0;k<gg;k++){
                if(*p1j==*p2j && *p1j!=-1){
                    repeat_value++;
                }else{
                    break;
                }
                if(k<gg-1){
                    p1j++;
                    p2j++;
                }
            }
            if(repeat_value==gg){
                repeat_count++;
                //printf("debug_1\n");
                *pr=j;
                //printf("debug_2\n");
                pr++;
                for(int* pcj=*p2i;pcj<*p2i+gg;pcj++){
                    *pcj=-1;
                }
            }


            if(j<ll-1){
                p2i++;
            }
        }

        if(i<ll-1){
            p1i++;
        }
    }

    int* repeat_record=(int*)malloc(repeat_count*sizeof(int));
    int* pi1=repeat_record;
    for(int* pi=temp_record;pi<temp_record+repeat_count;pi++,pi1++){
        *pi1=*pi;
    }
    free(temp_record);
    temp_record=NULL;

    if(repeat_count){
        //printf("repeat:%d\n",repeat_count);
    }

    *repeat_result=repeat_count;
    return repeat_record;
}


int function_i_find(int nn,int gg,double judge){
    //通过Pi的近似分布，递归进行计算（本处为单个函数，通过反复调用此函数实现），以随机数获得超连接中的每个个体

    double temp_i_middle=(double)(nn-gg+1)*(1.-pow(1.-judge,1./gg));
    int origin_i=(int)temp_i_middle+1-1;//进位和平移取值

    return origin_i;//输出范围[0,nn-1],程序所用数值为代数值减一
}


/// It is a bad way to name files with following functions, but it works, so I did not change it


char* function_file_name_txt_density(int xx){
    //使用整型数x转换成字符串“density\x.txt",仅限Windows使用（末尾反斜杠分隔符）
    int n10=0;
    int x1=xx;
    while(x1){
        x1=(int)x1/10;
        n10++;
    }

    char* c1=(char*)malloc((n10+1)*sizeof(char));
    //itoa(xx,c1,10);
    sprintf(c1,"%d",xx);

    char* c2=(char*)malloc((8+5+8)*sizeof(char));

    for(char* pi=c2+8;pi<c2+8+8-n10;pi++){
        *pi=48;
    }

    char* pj=c1;
    for(char* pi=c2+8+8-n10;pi<c2+8+8;pi++){
        *pi=*(pj++);
    }

    char* pi=c2;
    *(pi++)=100;
    *(pi++)=101;
    *(pi++)=110;
    *(pi++)=115;
    *(pi++)=105;
    *(pi++)=116;
    *(pi++)=121;
    *(pi++)=divide_sign;//"density\"

    pi=c2+8+8;
    *(pi++)=46;
    *(pi++)=116;
    *(pi++)=120;
    *(pi++)=116;
    *pi=0;//".txt\0"

    free(c1);
    c1=NULL;

    return c2;
}


char* function_file_name_txt_time_series(int group_n,int repeat_n){
    //为整洁起见并保证够用，组数与重复数均限制为4位十进制数

    group_n++;
    repeat_n++;
    //0开始变1开始

    char* c1=(char*)malloc(30*sizeof(char));

    int n10=0;
    int x1=group_n;
    while(x1){
        x1=(int)x1/10;
        n10++;
    }
    char* gg=(char*)malloc(5*sizeof(char));
    //itoa(group_n,gg,10);
    sprintf(gg,"%d",group_n);

    int n20=0;
    int x2=repeat_n;
    while(x2){
        x2=(int)x2/10;
        n20++;
    }
    char* rr=(char*)malloc(5*sizeof(char));
    //itoa(repeat_n,rr,10);
    sprintf(rr,"%d",repeat_n);

    char* pi=c1;
    *(pi++)=116;
    *(pi++)=105;
    *(pi++)=109;
    *(pi++)=101;
    *(pi++)=95;
    *(pi++)=115;
    *(pi++)=101;
    *(pi++)=114;
    *(pi++)=105;
    *(pi++)=101;
    *(pi++)=115;
    *(pi++)=divide_sign;
    //"time_series/"

    *(pi++)=103;//g
    for(int i=0;i<4-n10;i++){
        *(pi++)=48;
    }
    char* pm=gg;
    for(int i=0;i<n10;i++){
        *(pi++)=*(pm++);
    }

    *(pi++)=95;
    *(pi++)=114;//_r
    for(int i=0;i<4-n20;i++){
        *(pi++)=48;
    }
    pm=rr;
    for(int i=0;i<n20;i++){
        *(pi++)=*(pm++);
    }

    *(pi++)=46;
    *(pi++)=116;
    *(pi++)=120;
    *(pi++)=116;
    *pi=0;//".txt\0"

    free(gg);
    free(rr);


    return c1;
}


char* function_file_name_txt_hypergraph_(int group_n,int repeat_n){
    //为整洁起见并保证够用，组数与重复数均限制为4位十进制数

    group_n++;
    repeat_n++;
    //0开始变1开始

    char* c1=(char*)malloc(30*sizeof(char));

    int n10=0;
    int x1=group_n;
    while(x1){
        x1=(int)x1/10;
        n10++;
    }
    char* gg=(char*)malloc(5*sizeof(char));
    //itoa(group_n,gg,10);
    sprintf(gg,"%d",group_n);

    int n20=0;
    int x2=repeat_n;
    while(x2){
        x2=(int)x2/10;
        n20++;
    }
    char* rr=(char*)malloc(5*sizeof(char));
    //itoa(repeat_n,rr,10);
    sprintf(rr,"%d",repeat_n);

    char* pi=c1;
    *(pi++)=104;
    *(pi++)=121;
    *(pi++)=112;
    *(pi++)=101;
    *(pi++)=114;
    *(pi++)=103;
    *(pi++)=114;
    *(pi++)=97;
    *(pi++)=112;
    *(pi++)=104;
    *(pi++)=95;
    *(pi++)=divide_sign;
    //"hypergraph_/"

    *(pi++)=103;//g
    for(int i=0;i<4-n10;i++){
        *(pi++)=48;
    }
    char* pm=gg;
    for(int i=0;i<n10;i++){
        *(pi++)=*(pm++);
    }

    *(pi++)=95;
    *(pi++)=114;//_r
    for(int i=0;i<4-n20;i++){
        *(pi++)=48;
    }
    pm=rr;
    for(int i=0;i<n20;i++){
        *(pi++)=*(pm++);
    }

    *(pi++)=46;
    *(pi++)=116;
    *(pi++)=120;
    *(pi++)=116;
    *pi=0;//".txt\0"

    free(gg);
    free(rr);

    return c1;
}


int function_folder_create(){
    char* c1=(char*)malloc(8*sizeof(char));
    char* pi=c1;
    *(pi++)=100;
    *(pi++)=101;
    *(pi++)=110;
    *(pi++)=115;
    *(pi++)=105;
    *(pi++)=116;
    *(pi++)=121;
    *pi=0;
    mkdir(c1, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //创建输出结果的文件夹density
    //存储不同时刻分布，旧代码重复利用遗留问题——文件夹命名有错，暂未修改
    //需要开启迭代函数中的一个代码块

    char* c2=(char*)malloc(12*sizeof(char));
    pi=c2;
    *(pi++)=116;
    *(pi++)=105;
    *(pi++)=109;
    *(pi++)=101;
    *(pi++)=95;
    *(pi++)=115;
    *(pi++)=101;
    *(pi++)=114;
    *(pi++)=105;
    *(pi++)=101;
    *(pi++)=115;
    *pi=0;
    mkdir(c2, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //创建文件夹“time_series”用于存储所有组的时间序列文件

    char* c4=(char*)malloc(12*sizeof(char));
    pi=c4;
    *(pi++)=104;
    *(pi++)=121;
    *(pi++)=112;
    *(pi++)=101;
    *(pi++)=114;
    *(pi++)=103;
    *(pi++)=114;
    *(pi++)=97;
    *(pi++)=112;
    *(pi++)=104;
    *(pi++)=95;
    *pi=0;
    mkdir(c4, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    //创建文件夹“hypergraph_”用于存储所有组的超图结构

    free(c1);
    c1=NULL;
    free(c2);
    c2=NULL;
    free(c4);
    c4=NULL;

    return 0;
}
