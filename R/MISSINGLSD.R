#' Missing Plot in Latin Square Design(LSD)
#' @export
#' @param data a dataframe which contains values of LSD ,Index No.of row ,Index No.of column,Treatments/Index No.of Treatments in the 1st,2nd,3rd and 4th column respectively. In this dataframe, we will replace the missing value with 0 .
#' @param r the index no. of row containing the missing value.
#' @param c the index no. of column  containing the missing value.
#' @param t the index no. of Treatment containing the missing value.
#' @author  Arnab Roy , Debarghya Baul.
#' @importFrom stats aov qf
#' @description This function analyses LSD when there is one missing observation.
#' @details In LSD setup , if there is one missing observation we can use this function to estimate the missing observation along with Sum of Squares for testing the differential effect of the treatments. Here we estimate the missing observation twice by minimizing the SSE of the design.
#' @return x.hat : the least sqaure estimate of the missing observation.
#' @return  SSE.x.hat : Sum of Squares of Error of x.hat.
#' @return x.double.hat : the least square estimate of the missing observation under the null hypothesis , H0.
#' @return SSE.x.double.hat : Sum of Squares of Error of x.double.hat.
#' @return F.stat : Observed value of the Test Statistic.
#' @return F.crit.value : Critical value of the Test Statistic.
#' @examples  d=OrchardSprays
#' @examples d[8,1]=0
#' @examples  Missing.LSD(d,8,1,1)



Missing.LSD=function(data,r,c,t){
  if(is.character(data[,4])==TRUE){
    C_1=sum(data[data[,3]==c,1])
    R_1=sum(data[data[,2]==r,1])
    G=sum(data[,1])
    T_1=sum(data[data[,4]==LETTERS[t],1])
    v=length(unique(data[,4]))
    x.hat=(v*(R_1+C_1+T_1)-2*G)/((v-1)*(v-2))
    data[data[,2]==r&data[,3]==c&data[,4]==LETTERS[t],1]=x.hat

    y=data[,1]
    row=as.factor(data[,2])
    col=as.factor(data[,3])
    trt=as.factor(data[,4])
    s=summary(aov(y~row+col+trt))
    s1=s[[1]][2]$`Sum Sq`
    SSE_x=s1[4]

    x.double.hat=(v*(R_1+C_1)-G)/(v-1)^2;
    data[data[,2]==r&data[,3]==c&data[,4]==LETTERS[t],1]=x.double.hat

    y1=data[,1]
    s3=summary(aov(y1~row+col))
    s2=s3[[1]][2]$`Sum Sq`
    SSE.x.double.hat=s2[3]

    F_stat=((SSE.x.double.hat-SSE_x)/(v-1))/(SSE_x/((v-1)*(v-2)-1))
    F.crit=qf(0.95,v-1,(v-1)*(v-2)-1)

    list(x.hat=x.hat,SSE.x.hat=SSE_x,x.double.hat=x.double.hat,SSE.x.double.hat=
           SSE.x.double.hat,F.stat=F_stat,F.crit.value=F.crit)
  }else
  {
    C_1=sum(data[data[,3]==c,1])
    R_1=sum(data[data[,2]==r,1])
    G=sum(data[,1])
    T_1=sum(data[data[,4]==t,1])
    v=length(unique(data[,4]))
    x.hat=(v*(R_1+C_1+T_1)-2*G)/((v-1)*(v-2))
    data[data[,2]==r&data[,3]==c&data[,4]==t,1]=x.hat

    y=data[,1]
    row=as.factor(data[,2])
    col=as.factor(data[,3])
    trt=as.factor(data[,4])
    s=summary(aov(y~row+col+trt))
    s1=s[[1]][2]$`Sum Sq`
    SSE_x=s1[4]

    x.double.hat=(v*(R_1+C_1)-G)/(v-1)^2;
    data[data[,2]==r&data[,3]==c&data[,4]==t,1]=x.double.hat

    y1=data[,1]
    s3=summary(aov(y1~row+col))
    s2=s3[[1]][2]$`Sum Sq`
    SSE.x.double.hat=s2[3]

    F_stat=((SSE.x.double.hat-SSE_x)/(v-1))/(SSE_x/((v-1)*(v-2)-1))
    F.crit=qf(0.95,v-1,(v-1)*(v-2)-1)

    list(x.hat=x.hat,SSE.x.hat=SSE_x,x.double.hat=x.double.hat,SSE.x.double.hat=
           SSE.x.double.hat,F.stat=F_stat,F.crit.value=F.crit)
  }

}
