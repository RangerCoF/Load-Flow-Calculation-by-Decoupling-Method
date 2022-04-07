/*
基于Eigen的潮流计算程序
潮流计算方法采用PQ分解法
模型使用IEEE14结点测试系统 数据来源: https://labs.ece.uw.edu/pstca/pf14/pg_tca14bus.htm
编写时考虑了代码灵活性，如果想要改变模型参数和结点关系，可在注释有符号"%"处修改
Designed by Zhixuan Ge, CAU
*/

//通过标准库std::complex实现复数运算
//通过开源矩阵运算库Eigen实现矩阵运算
#include <iostream>         //标准输入输出
#include<cstring>           //字符串
#include<complex>           //标准复数
#include"eigen/Eigen/Dense"       //矩阵运算

#define complexd complex<double>    //对复数类进行宏定义方便编写
#define epsilon 1e-6                //收敛判定条件设为10^(-6)

using namespace std;                //设定命名空间
using namespace Eigen;

void Shrink(VectorXd& V, VectorXi& Vseq);               //向量维度缩减函数
void Restore(VectorXd& V, VectorXi& Vseq, int n);       //向量维度恢复函数

int main()
{
    int NumNode = 14;//输入结点数    %

    //建立并初始化结点导纳矩阵  %
    MatrixXd G = MatrixXd::Constant(NumNode, NumNode, 0), B = MatrixXd::Constant(NumNode, NumNode, 0);
    G(0, 0) = 5.94327, G(0, 1) = G(1, 0) = -4.99913, G(0, 4) = G(4, 0) = -0.94414;
    G(1, 1) = 9.52132, G(1, 2) = G(2, 1) = -1.13502, G(1, 3) = G(3, 1) = -1.68603, G(1, 4) = G(4, 1) = -1.70114;
    G(2, 2) = 3.12099, G(2, 3) = G(3, 2) = -1.98598;
    G(3, 3) = 10.51299, G(3, 4) = G(4, 3) = -6.84098;
    G(4, 4) = 9.48626;
    G(5, 5) = 8.1748, G(5, 10) = G(10, 5) = -1.95503, G(5, 11) = G(11, 5) = -3.12084, G(5, 12) = G(12, 5) = -3.09893;
    G(8, 8) = 5.32606, G(8, 9) = G(9, 8) = -3.90205, G(8, 13) = G(13, 8) = -1.42401;
    G(9, 9) = 5.78293, G(9, 10) = G(10, 9) = -1.88088;
    G(10, 10) = 3.83591;
    G(11, 11) = 5.60986, G(11, 12) = G(12, 11) = -2.48902;
    G(12, 12) = 6.72495, G(12, 13) = G(13, 12) = -1.13699;
    G(13, 13) = 2.561;

    B(0, 0) = -19.2843, B(0, 1) = B(1, 0) = 15.26309, B(0, 4) = B(4, 0) = 4.07221;
    B(1, 1) = -30.27072, B(1, 2) = B(2, 1) = 4.78186, B(1, 3) = B(3, 1) = 5.11584, B(1, 4) = B(4, 1) = 5.19393;
    B(2, 2) = -9.81148, B(2, 3) = B(3, 2) = 5.06882;
    B(3, 3) = -38.63517, B(3, 4) = B(4, 3) = 21.57855, B(3, 6) = B(6, 3) = 4.88951, B(3, 8) = B(8, 3) = 1.8555;
    B(4, 4) = -35.36477, B(4, 5) = B(5, 4) = 4.25754;
    B(5, 5) = -18.12098, B(5, 10) = B(10, 5) = 4.09407, B(5, 11) = B(11, 5) = 3.95651, B(5, 12) = B(12, 5) = 6.10276;
    B(6, 6) = -19.549, B(6, 7) = B(7, 6) = 5.67698, B(6, 8) = B(8, 6) = 9.09008;
    B(7, 7) = -5.67698;
    B(8, 8) = -24.09251, B(8, 9) = B(9, 8) = 10.36539, B(8, 13) = B(13, 8) = 3.02905;
    B(9, 9) = -14.76834, B(9, 10) = B(10, 9) = 4.40294;
    B(10, 10) = -8.49702;
    B(11, 11) = -6.20819, B(11, 12) = B(12, 11) = 2.25197;
    B(12, 12) = -10.66969, B(12, 13) = B(13, 12) = 2.31496;
    B(13, 13) = -5.34401;

    MatrixXcd Y = G + complexd(0, 1) * B;
    //
    
    //记录各结点数量并初始化结点列表   %
    int  NumNodePV = 4, NumNodePQ = 9;
    int  NumNodeNotR = NumNode - 1, NumNodeR = 1;
    VectorXi NodeR(NumNodeR), NodePV(NumNodePV), NodePQ(NumNodePQ), NodeNotR(NumNode - 1);
    NodeR << 0;
    NodePV << 1, 2, 5, 7;
    NodePQ << 3, 4, 6, 8, 9, 10, 11, 12, 13;
    NodeNotR << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13;
    //
    
    //记录各结点电压并初始化   %
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
            cout << "迭代次数: " << i + 1 << endl << "精度: " << epsilon << "\n" << endl;    //输出迭代次数
            break;
        }
        //
    }

    //输出结果
    cout << "U: " << endl << U << "\n" << endl;                           //输出U
    cout << "delta: " << endl << delta / 3.1415926 * 180 << "\n" << endl;     //以角度形式输出δ和δ

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