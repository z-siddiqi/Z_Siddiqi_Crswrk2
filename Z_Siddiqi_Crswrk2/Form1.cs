using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Z_Siddiqi_Crswrk2
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            //Initialise the text boxes
            textBox1.Text = "9";
            textBox2.Text = "1.5";
            textBox3.Text = "1";
            textBox4.Text = "5.93";
            textBox5.Text = "0";
            textBox6.Text = "0";
            textBox7.Text = "4";
        }

        // Calculate button
        private void button1_Click(object sender, EventArgs e)
        {
            double b = Convert.ToDouble(textBox1.Text);
            double cr = Convert.ToDouble(textBox2.Text);
            double taperRatio = Convert.ToDouble(textBox3.Text);
            double ae = Convert.ToDouble(textBox4.Text);
            double alphaZeroLift = Convert.ToDouble(textBox5.Text) * Math.PI / 180;
            double washout = Convert.ToDouble(textBox6.Text);

            double aoaRoot = Convert.ToDouble(textBox7.Text);

            double[] RHS = new double[6];
            double[][] LHS = MatrixCreate(6, 6);

            for (int i = 1; i <= 6; i++)
            {
                int column = 0;

                double alpha = (aoaRoot - washout / 7 * i) * Math.PI / 180;

                double phi = Math.Acos((double)i / 7);
                double mu = Mu(ae, b, cr, taperRatio, phi);

                RHS[i - 1] = mu * (alpha - alphaZeroLift) * Math.Sin(phi);

                for (int j = 1; j <= 11; j += 2)
                {
                    LHS[i - 1][column] = Math.Sin(j * phi) * (j * mu + Math.Sin(phi));

                    column++;
                }
            }

            double[][] ILHS = MatrixInverse(LHS);
            double[] A = MatrixVectorProduct(ILHS, RHS);

            double CL = A[0] * Math.PI * Ar(b, cr, taperRatio);
            double CD = Math.Pow(CL, 2) / (Math.PI * Ar(b, cr, taperRatio)) *
                    (1 + 3 * (Math.Pow(A[1], 2) / Math.Pow(A[0], 2)) + 
                    5 * (Math.Pow(A[2], 2) / Math.Pow(A[0], 2)) +
                    7 * (Math.Pow(A[3], 2) / Math.Pow(A[0], 2)) +
                    9 * (Math.Pow(A[4], 2) / Math.Pow(A[0], 2)) +
                    11 * (Math.Pow(A[5], 2) / Math.Pow(A[0], 2)));

            // Set the relevant text boxes to their calculated value
            textBox8.Text = Convert.ToString(A[0]);
            textBox9.Text = Convert.ToString(A[1]);
            textBox10.Text = Convert.ToString(A[2]);
            textBox11.Text = Convert.ToString(A[3]);
            textBox12.Text = Convert.ToString(CL);
            textBox13.Text = Convert.ToString(CD);
        }

        // Exit button
        private void button2_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        static double Mu(double ae, double b, double cr, double taperRatio, double phi)
        {
            double AR = Ar(b, cr, taperRatio);

            return (ae / (2 * AR * (1 + taperRatio))) * (1 + Math.Cos(phi) * (taperRatio - 1));
        }

        static double Ar(double b, double cr, double taperRatio)
        {
            return 2 * b / (cr * (1 + taperRatio));
        }

        // #############################################################//
        //                 matrix methods                               //
        //##############################################################//

        static double[][] MatrixCreate(int rows, int cols)
        {
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols];
            return result;
        }

        static double[][] MatrixDuplicate(double[][] matrix)
        {
            // allocates/creates a duplicate of a matrix.
            double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; ++i) // copy the values
                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            return result;
        }

        static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
        {
            // Doolittle LUP decomposition with partial pivoting.
            // rerturns: result is L (with 1s on diagonal) and U;
            // perm holds row permutations; toggle is +1 or -1 (even or odd)
            int rows = matrix.Length;
            int cols = matrix[0].Length; // assume square
            if (rows != cols)
                throw new Exception("Attempt to decompose a non-square m");

            int n = rows; // convenience

            double[][] result = MatrixDuplicate(matrix);

            perm = new int[n]; // set up row permutation result
            for (int i = 0; i < n; ++i) { perm[i] = i; }

            toggle = 1; // toggle tracks row swaps.
                        // +1 -greater-than even, -1 -greater-than odd. used by MatrixDeterminant

            for (int j = 0; j < n - 1; ++j) // each column
            {
                double colMax = Math.Abs(result[j][j]); // find largest val in col
                int pRow = j;

                // reader Matt V needed this:
                for (int i = j + 1; i < n; ++i)
                {
                    if (Math.Abs(result[i][j]) > colMax)
                    {
                        colMax = Math.Abs(result[i][j]);
                        pRow = i;
                    }
                }
                // Not sure if this approach is needed always, or not.

                if (pRow != j) // if largest value not on pivot, swap rows
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[pRow]; // and swap perm info
                    perm[pRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }

                // --------------------------------------------------
                // This part added later (not in original)
                // and replaces the 'return null' below.
                // if there is a 0 on the diagonal, find a good row
                // from i = j+1 down that doesn't have
                // a 0 in column j, and swap that good row with row j
                // --------------------------------------------------

                if (result[j][j] == 0.0)
                {
                    // find a good row to swap
                    int goodRow = -1;
                    for (int row = j + 1; row < n; ++row)
                    {
                        if (result[row][j] != 0.0)
                            goodRow = row;
                    }

                    if (goodRow == -1)
                        throw new Exception("Cannot use Doolittle's method");

                    // swap rows so 0.0 no longer on diagonal
                    double[] rowPtr = result[goodRow];
                    result[goodRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[goodRow]; // and swap perm info
                    perm[goodRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }

                for (int i = j + 1; i < n; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < n; ++k)
                    {
                        result[i][k] -= result[i][j] * result[j][k];
                    }
                }
            } // main j column loop
            return result;
        }// MatrixDecompose

        static double[][] MatrixInverse(double[][] matrix)
        {
            int n = matrix.Length;
            double[][] result = MatrixDuplicate(matrix);

            int[] perm;
            int toggle;
            double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
            if (lum == null)
                throw new Exception("Unable to compute inverse");

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }

                double[] x = HelperSolve(lum, b); // 

                for (int j = 0; j < n; ++j)
                    result[j][i] = x[j];
            }
            return result;
        }

        static double[] MatrixVectorProduct(double[][] matrix, double[] vector)
        {
            // result of multiplying an n x m matrix by a m x 1 
            // column vector (yielding an n x 1 column vector)
            int mRows = matrix.Length; int mCols = matrix[0].Length;
            int vRows = vector.Length;
            if (mCols != vRows)
                throw new Exception("Non-conformable matrix and vector");
            double[] result = new double[mRows];
            for (int i = 0; i < mRows; ++i)
                for (int j = 0; j < mCols; ++j)
                    result[i] += matrix[i][j] * vector[j];
            return result;
        }

        static double[] HelperSolve(double[][] luMatrix, double[] b)
        {
            // before calling this helper, permute b using the perm array
            // from MatrixDecompose that generated luMatrix
            int n = luMatrix.Length;
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum / luMatrix[i][i];
            }
            return x;
        }
    }
}
