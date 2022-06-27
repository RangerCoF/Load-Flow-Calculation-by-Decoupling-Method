/*
基于Eigen的潮流计算程序
潮流计算方法采用PQ分解法
模型使用IEEE14结点测试系统 数据来源: https://labs.ece.uw.edu/pstca/pf14/pg_tca14bus.htm
编写时考虑了代码灵活性，如果想要改变模型参数和结点关系，可在注释有符号"%"处修改
Designed by Zhixuan Ge, CAU
Github: https://github.com/CallofFood/Load-Flow-Calculation-by-Decoupling-Method
*/

//通过标准库std::complex实现复数运算
//通过开源矩阵运算库Eigen实现矩阵运算
#include <iostream>         //标准输入输出
#include<cstring>           //字符串
#include<complex>           //标准复数
#include"eigen/Eigen/Dense"       //矩阵运算

#define complexd complex<double>    //对复数类进行宏定义方便编写
#define epsilon 1e-6                //收敛判定条件设为10^(-6)


using namespace Eigen;
using namespace std;                //设定命名空间

void Shrink(VectorXd& V, VectorXi& Vseq);               //向量维度缩减函数
void Restore(VectorXd& V, VectorXi& Vseq, int n);       //向量维度恢复函数
void CmpOut(complexd c);                                //格式化输出复数

int main()
{
    int NumNode = 14;//输入结点数    %
    int NumLine = 20;//路线数        %

    //初始化路线参数 %
    VectorXi linehead(NumLine), linetail(NumLine);              //线路两端对应结点
    VectorXd Gline(NumLine), Ghead(NumLine), Gtail(NumLine);    //线路导纳和两端对地导纳
    VectorXd Bline(NumLine), Bhead(NumLine), Btail(NumLine);
    VectorXd GnodeN(14), BnodeN(14);                            //结点对地导纳

    linehead << 1, 2, 2, 1, 2, 3, 4, 7, 7, 9, 6, 6, 6, 9, 10, 12, 13, 5, 4, 4;
    linetail << 2, 3, 4, 5, 5, 4, 5, 8, 9, 10, 11, 12, 13, 14, 11, 13, 14, 6, 7, 9;
    linehead = linehead.array() - 1;
    linetail = linetail.array() - 1;

    Gline << 4.99913, 1.13502, 1.68603, 0.94414, 1.70114, 1.98598, 6.84098, 0.00000, 0.00000, 3.90205, 1.95503, 3.12084, 3.09893, 1.42401, 1.88088, 2.48902, 1.13699, 0.00000, 0.00000, 0.00000;
    Ghead = Gtail = VectorXd::Constant(20, 0);

    Bline << -15.26309, -4.78186, -5.11584, -4.07221, -5.19393, -5.06882, -21.57855, -5.67698, -9.09008, -10.36539, -4.09407, -3.95621, -6.10276, -3.02905, -4.40294, -2.25197, -2.31496, -4.25745, -4.88951, -1.85550;
    Bhead = Btail = VectorXd::Constant(20, 0);
    Bhead.head(7) << 0.0264, 0.0219, 0.0187, 0.0246, 0.017, 0.0173, 0.0064;
    Bhead.tail(3) << -0.31063, -0.10999, -0.05936;
    Btail.head(7) << 0.0264, 0.0219, 0.0187, 0.0246, 0.017, 0.0173, 0.0064;
    Btail.tail(3) << 0.28951, 0.10757, 0.05752;

    GnodeN = BnodeN = VectorXd::Constant(NumNode, 0);
    BnodeN(8) = 0.19;
    //
    
    //初始化结点分类参数   %
    int  NumNodePV = 4, NumNodePQ = 9;
    int  NumNodeNotR = NumNode - 1, NumNodeR = 1;
    VectorXi NodeR(NumNodeR), NodePV(NumNodePV), NodePQ(NumNodePQ), NodeNotR(NumNode - 1);
    NodeR << 0;
    NodePV << 1, 2, 5, 7;
    NodePQ << 3, 4, 6, 8, 9, 10, 11, 12, 13;
    NodeNotR << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13;
    //
    
    //初始化结点电压参数   %
    VectorXd P(NumNode), Q(NumNode);
    VectorXd  U = VectorXd::Constant(NumNode, 1), delta = VectorXd::Constant(NumNode, 0);
    P << 2.324, 0.183, -0.942, -0.478, -0.076, -0.112, 0, 0, -0.295, -0.09, -0.035, -0.061, -0.135, -0.149;
    Q << -0.169, 0.297, 0.044, 0.039, -0.016, 0.047, 0, 0.174, -0.166, -0.058, -0.018, -0.016, -0.058, -0.05;
    U(0) = 1.06;
    U(1) = 1.045;
    U(2) = 1.01;
    U(5) = 1.07;
    U(7) = 1.09;
    //

    //由路线参数生成结点电压矩阵
    MatrixXcd Yline = Gline + complexd(0, 1) * Bline;
    MatrixXcd Yhead = Ghead + complexd(0, 1) * Bhead;
    MatrixXcd Ytail = Gtail + complexd(0, 1) * Btail;
    MatrixXcd YnodeN = GnodeN + complexd(0, 1) * BnodeN;

    MatrixXd  B(NumNode, NumNode), G(NumNode, NumNode);
    MatrixXcd Y(NumNode, NumNode);
    Y = B = G = MatrixXd::Constant(NumNode, NumNode, 0);


    for (int i = 0; i < NumLine; i++)
    {
        Y(linehead(i), linetail(i)) = Y(linetail(i), linehead(i)) = -Yline(i);
        Y(linehead(i), linehead(i)) += (Yline(i) + Yhead(i));
        Y(linetail(i), linetail(i)) += (Yline(i) + Ytail(i));
    }

    for (int i = 0; i < NumNode; i++)
    {
        Y(i, i) += YnodeN(i);
    }

    G = Y.array().real();
    B = Y.array().imag();
    //
    

    //获取B'和B''
    //B'为B删去和参考节点对应的行和列，B''为B仅保留和PQ结点相关的行和列
    MatrixXd B1(NumNode - 1, NumNode - 1), B2(NumNodePQ, NumNodePQ);
    MatrixXd tempQ;

    tempQ = MatrixXd::Constant(NumNode - 1, NumNode, 0);
    for (int i = 0; i < NumNode - 1; i++)
    {
        tempQ.row(i) = B.row(NodeNotR(i));
    }
    for (int i = 0; i < NumNode - 1; i++)
    {
        B1.col(i) = tempQ.col(NodeNotR(i));
    }

    tempQ = MatrixXd::Constant(NumNodePQ, NumNode, 0);
    for (int i = 0; i < NumNodePQ; i++)
    {
        tempQ.row(i) = B.row(NodePQ(i));
    }
    for (int i = 0; i < NumNodePQ; i++) 
    {
        B2.col(i) = tempQ.col(NodePQ(i));
    }

    B1 = -B1.inverse();     //每次迭代矩阵一致
    B2 = -B2.inverse();     //为提高效率直接求解逆矩阵
    //
    

    //建立求解相关矩阵
    VectorXd dP, dQ, ddelta, dU;
    VectorXd dP_U, dQ_U, Ubyddelta;
    VectorXd tempV = VectorXd::Constant(NumNode, 0);
    MatrixXd delta2d(NumNode,NumNode);
    MatrixXd Cosdelta, Sindelta;
    //

    //由于三角函数矩阵tri(i,j)中每个元素每次迭代都需要使用两次（该次迭代电压运算和下次迭代相角运算）
    //为减少三角函数运算次数，使用矩阵Cosδ和Sinδ储存对应的值
    //由初值计算第1次迭代时的三角函数值
    for (int i = 0; i < NumNode; i++)           //由向量delta(i)得到矩阵delta2d(i,j)
        for (int j = 0; j < NumNode; j++)
        {
            delta2d(i, j) = delta(i) - delta(j);
        }
    Cosdelta = delta2d.array().cos();           //三角矩阵运算
    Sindelta = delta2d.array().sin();
    //

    //迭代开始
    for (int i = 0; i < 200; i++)
    {
        //利用矩阵运算（矩阵乘法和同矩阵元素相乘）进行高效计算
        //理论上对于NumNode个结点，仅需要参考节点之外的NumNode-1个结点的dP进行迭代
        //但算法为了维度一致先计算出含参考节点的NumNode维dP，之后操作删去参考节点对应元素
        dP = P.array() - ((G.array() * Cosdelta.array() + B.array() * Sindelta.array()).matrix() * U).array() * U.array();
        dP_U = dP.array()/U.array();
        Shrink(dP_U,NodeNotR);          //删去参考节点对应元素（保留所有非参考结点元素）

        Ubyddelta = B1 * dP_U;                  //迭代                  
        Restore(Ubyddelta, NodeNotR, NumNode);  //恢复至NumNode维
        ddelta = Ubyddelta.array() / U.array(); //计算dδ
        delta += ddelta;                        //更新dδ

        //更新Cosδ和Sinδ
        for (int i = 0; i < NumNode; i++)
            for (int j = 0; j < NumNode; j++)
            {
                delta2d(i, j) = delta(i) - delta(j);
            }
        Cosdelta = delta2d.array().cos();
        Sindelta = delta2d.array().sin();
        //

        
        //同上，利用矩阵求解dQ
        dQ = Q.array() - ((G.array() * Sindelta.array() - B.array() * Cosdelta.array()).matrix() * U).array() * U.array();
        dQ_U = dQ.array() / U.array();
        
        Shrink(dQ_U, NodePQ);       //保留PQ结点相关元素，其他删除

        dU = B2 * dQ_U;             //迭代
        
        Restore(dU, NodePQ, NumNode);   //重建至NumNode维

        U += dU;                    //更新U
        
        //判断收敛
        Shrink(dP, NodeNotR);   //只需要非参考节点dP和PQ结点的dQ，其他dP和dQ没有意义且不参与迭代，且可能影响收敛判断
        Shrink(dQ, NodePQ);     //因此收敛判断前删去相关结点对应元素
        if (dP.array().abs().maxCoeff() < epsilon && dQ.array().abs().maxCoeff() < epsilon && dU.array().abs().maxCoeff() < epsilon)
        {
            cout << "Iteration Ordinal Number: " << i + 1 << endl << "Accuracy: " << epsilon << "\n" << endl;    //输出迭代次数
            break;
        }
        //
    }

    //输出迭代结果
    cout << "U: " << endl << U << "\n" << endl;                           //输出U
    cout << "delta: " << endl << delta / 3.1415926 * 180 << "\n" << endl;     //以角度形式输出δ和δ
    //

    //计算复功率
    VectorXcd Uph = U.array() * delta.array().cos() + complexd(0, 1) * U.array() * delta.array().sin();
    MatrixXcd S(NumLine,3);
    for (int i = 0; i < NumLine; i++)
    {
        S(i, 0) = Uph(linehead(i)) * conj(Uph(linehead(i)) * Yhead(i) + (Uph(linehead(i)) - Uph(linetail(i))) * Yline(i));
        S(i, 1) = Uph(linetail(i)) * conj(Uph(linetail(i)) * Ytail(i) + (Uph(linetail(i)) - Uph(linehead(i))) * Yline(i));
    }
    S.col(2) = S.col(0) + S.col(1);

    //输出复功率
    for (int i = 0; i < NumLine; i++)
    {
        cout << "head: " << linehead(i) + 1 << ", tail: " << linetail(i) + 1 << ", S = ";
        CmpOut(S(i, 0));
        cout << endl;
        
        cout << "head: " << linetail(i) + 1 << ", tail: " << linehead(i) + 1 << ", S = ";
        CmpOut(S(i, 1));
        cout << endl;

        cout << "head: " << linehead(i) + 1 << ", tail: " << linetail(i) + 1 << ", dS = ";
        CmpOut(S(i, 2));
        cout << endl;
        cout << "\n";
    }
    //cout << S;

    return 0;
}

//该函数用于删去向量V中与Vseq储存的结点无关的元素
void Shrink(VectorXd & V, VectorXi & Vseq)
{
    int size = Vseq.size();
    VectorXd tempV = VectorXd::Constant(size, 0);
    for (int i = 0; i < size; i++)
        tempV(i) = V(Vseq(i));
    V = tempV;
}

//该函数用于恢复向量V恢复至n维，原本V中元素对应位置由Vseq决定
void Restore(VectorXd& V, VectorXi& Vseq, int n)
{
    int size = Vseq.size();
    VectorXd tempV = VectorXd::Constant(n, 0);
    for (int i = 0; i < size; i++)
        tempV(Vseq(i)) = V(i);
    V = tempV;
}

void CmpOut(complexd c)
{
    string PorN;
    double cImag;
    if (c.imag() > 0)
    {
        PorN = " + ";
        cImag = c.imag();
    }
    else
    {
        PorN = " - ";
        cImag = -c.imag();
    }
    cout << c.real() << PorN << "j" << cImag;
}