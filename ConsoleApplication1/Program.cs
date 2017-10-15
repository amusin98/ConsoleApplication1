using System;

namespace ConsoleApplication1
{
    class Program
    {
        static void Main(string[] args)
        {
            // TM - коэффициент теплопроводности 
            double TM = 100;
            double r = 1;  //радиус сечения
            double h = 10;  //быстрота теплообмена (конвекция)
            double q = 5; // плотность теплового потока
            double tOut = 40; //температура окружающей среды
            int kEls = 5;  //количество кон. элементов
            int Ui = 5; // количество строк
            double dt = 1;  //шаг по времени
            // kNodes - количество узлов
            int kNodes = kEls+1;
            // lEl - длина одного элемента
            double lEl = 10;
            // sCut - размер сечения (площадь поперечного сечения)
            double sCut = Math.PI * Math.Pow(r, 2);
            double CPM = 15; //теплоемк * плотн 
          

            
            double[,] T = new double[Ui, kNodes]; //вектор узловых значений температуры (решение)

            double[,] F = new double[Ui, kNodes]; //вектор нагрузки
            double[] midF = new double[kNodes];//вектор Fcp
            double[,] B = new double[kNodes, kNodes];
            double[,] P = new double[kNodes, kNodes];
            double[,] K = new double[kNodes, kNodes]; // матрица теплопроводности(глобальная)
            double[,] C = new double[kNodes, kNodes]; // матрица демпфирования(глобальная)


            int left = 100; //левые граничные условия
            int right = 500; //правые граничные условия

            T[0, 0] = left;
            T[0, kNodes-1] = right;

            double b = left;
            double a = (right - b) / ((kNodes-1)*lEl);
            //Шаг по х
            double xx = lEl;

            // Начальные условия
            for (int i = 1; i < kNodes - 1; i++)
            {
                if (i == 1)
                {
                    T[0, i] = a * xx + b;
                }
                else
                {
                    xx = xx + lEl;
                    T[0, i] = a * xx + b;
                }

            }

            // Вектор нагрузки для нулевого момента времени
            F[0, 0] = F[0, 1] = -(h * (left - tOut) * sCut) / 2;
            F[0, F.GetLength(1) - 1] = F[0, F.GetLength(1) - 2] = q * sCut / 2;

            // Локальные матрицы K и C
            double[,] kTTmp = new double[2, 2];
            double[,] kATmp = new double[2, 2];
            double[,] cTTmp = new double[2, 2];
            double[,] cATmp = new double[2, 2];
            kTTmp[0, 0] = sCut * TM / lEl;
            kTTmp[0, 1] = -kTTmp[0, 0];
            kTTmp[1, 0] = -kTTmp[0, 0];
            kTTmp[1, 1] = kTTmp[0, 0];
            cTTmp[0, 0] = CPM * sCut * lEl / 3;
            cTTmp[0, 1] = cTTmp[0, 0] / 2;
            cTTmp[1, 0] = cTTmp[0, 0] / 2;
            cTTmp[1, 1] = cTTmp[0, 0];

            // Заполнение глобальных матриц 
            for (int i = 0; i < kNodes - 1; i++)
            {
                
                    K[i, i] += kTTmp[0, 0];
                    K[i, i + 1] += kTTmp[0, 1];
                    K[i + 1, i] += kTTmp[1, 0];
                    K[i + 1, i + 1] += kTTmp[1, 1];
                    C[i, i] += cTTmp[0, 0];
                    C[i, i + 1] += cTTmp[0, 1];
                    C[i + 1, i] += cTTmp[1, 0];
                    C[i + 1, i + 1] += cTTmp[1, 1];
            }

          
             for (int i = 0; i < kNodes; i++)
                for (int j = 0; j < kNodes; j++)
                    C[i, j] *= 2 / dt;

    

            // Формирование матриц для уравнения
            for (int i = 0; i < kNodes; i++)
                for (int j = 0; j < kNodes; j++)
                {
                    B[i, j] = K[i, j] + C[i, j];
                    P[i, j] = C[i, j] - K[i, j];
                }

            // Заполнение матрицы температур
            for (int tNumb = 1; tNumb < T.GetLength(0); tNumb++)
            {
                // Вектор нагрузки
                F[tNumb, 0] = -(h * (T[tNumb - 1, 0] - tOut) * sCut) / 2;
                F[tNumb, 1] = -(h * (T[tNumb - 1, 1] - tOut) * sCut) / 2;
                F[tNumb, F.GetLength(1) - 1] = F[tNumb, F.GetLength(1) - 2] = q * sCut / 2;

                for (int j = 0; j < kNodes; j++)
                {
                    midF[j] = F[tNumb - 1, j] + F[tNumb, j];
                }

                //Матрица результат умеожения матрицы Р на вектор Т
                double[,] tmpBVect = new double[kNodes, 1];
                //Вектор, содержащий температуры на предыдущем слое
                double[,] tmpT = new double[kNodes, 1];
                //Вектор решние СЛАУ
                double[] solve = new double[kNodes];
                //Столбец свободных членов для системы
                double[] bVect = new double[kNodes];

                for (int i = 0; i < kNodes; i++)
                    tmpT[i, 0] = T[tNumb - 1, i];

                tmpBVect = mult_matrix(P, tmpT);
                
                for (int i = 0; i < kNodes; i++)
                    bVect[i] = tmpBVect[i, 0] - midF[i];
                
                


                // Решение СЛАУ
                GausMethod method = new GausMethod((uint)kNodes, (uint)kNodes);
                //заполняем правую часть
                method.RightPart = bVect;

                //заполняем матрицу
                for (int i = 0; i < kNodes; i++)
                    for (int j = 0; j < kNodes; j++)
                        method.Matrix[i][j] = B[i, j];

                //решаем матрицу
                method.SolveMatrix();

                //сохраняем ответ

                for (int i = 0; i < kNodes; i++)
                {
                    solve[i] = method.Answer[i];
                }

                // Заносим значения для найденного момента времени в соответствующую строку матрицы
                for (int j = 0; j < kNodes; j++)
                    T[tNumb, j] = solve[j];
            }

            // Вывод результатов
                for (int i = 0; i < T.GetLength(0); i++)
                {
                    for (int j = 0; j < T.GetLength(1); j++)
                    {
                        Console.Write("{0:f2} ", T[i, j]);
                    }
                    Console.WriteLine();
                }

            Console.ReadKey();
        }
       

        // Умножение матриц
        public static double[,] mult_matrix(double[,] array1, double[,] arrary2)
        {
            int commonLength = array1.GetLength(1);

            double[,] res = new double[array1.GetLength(0), arrary2.GetLength(1)];

            for (int i = 0; i < res.GetLength(0); i++)
            {
                for (int j = 0; j < res.GetLength(1); j++)
                {
                    double nextVal = 0;
                    for (int k = 0; k < commonLength; k++)
                    {
                        nextVal += array1[i, k] * arrary2[k, j];
                    }
                    res[i, j] = nextVal;
                }
            }
            return res;
        }
    }


    class GausMethod
    {
        public uint RowCount;
        public uint ColumCount;
        public double[][] Matrix { get; set; }
        public double[] RightPart { get; set; }
        public double[] Answer { get; set; }

        public GausMethod(uint Row, uint Colum)
        {
            RightPart = new double[Row];
            Answer = new double[Row];
            Matrix = new double[Row][];
            for (int i = 0; i < Row; i++)
                Matrix[i] = new double[Colum];
            RowCount = Row;
            ColumCount = Colum;

            //обнулим массив
            for (int i = 0; i < Row; i++)
            {
                Answer[i] = 0;
                RightPart[i] = 0;
                for (int j = 0; j < Colum; j++)
                    Matrix[i][j] = 0;
            }
        }

        private void SortRows(int SortIndex)
        {

            double MaxElement = Matrix[SortIndex][SortIndex];
            int MaxElementIndex = SortIndex;
            for (int i = SortIndex + 1; i < RowCount; i++)
            {
                if (Matrix[i][SortIndex] > MaxElement)
                {
                    MaxElement = Matrix[i][SortIndex];
                    MaxElementIndex = i;
                }
            }

            //теперь найден максимальный элемент ставим его на верхнее место
            if (MaxElementIndex > SortIndex)//если это не первый элемент
            {
                double Temp;

                Temp = RightPart[MaxElementIndex];
                RightPart[MaxElementIndex] = RightPart[SortIndex];
                RightPart[SortIndex] = Temp;

                for (int i = 0; i < ColumCount; i++)
                {
                    Temp = Matrix[MaxElementIndex][i];
                    Matrix[MaxElementIndex][i] = Matrix[SortIndex][i];
                    Matrix[SortIndex][i] = Temp;
                }
            }
        }

        public int SolveMatrix()
        {
            if (RowCount != ColumCount)
                return 1; //нет решения

            for (int i = 0; i < RowCount - 1; i++)
            {
                SortRows(i);
                for (int j = i + 1; j < RowCount; j++)
                {
                    if (Matrix[i][i] != 0) //если главный элемент не 0, то производим вычисления
                    {
                        double MultElement = Matrix[j][i] / Matrix[i][i];
                        for (int k = i; k < ColumCount; k++)
                            Matrix[j][k] -= Matrix[i][k] * MultElement;
                        RightPart[j] -= RightPart[i] * MultElement;
                    }
                    //для нулевого главного элемента просто пропускаем данный шаг
                }
            }

            //ищем решение
            for (int i = (int)(RowCount - 1); i >= 0; i--)
            {
                Answer[i] = RightPart[i];

                for (int j = (int)(RowCount - 1); j > i; j--)
                    Answer[i] -= Matrix[i][j] * Answer[j];

                if (Matrix[i][i] == 0)
                    if (RightPart[i] == 0)
                        return 2; //множество решений
                    else
                        return 1; //нет решения

                Answer[i] /= Matrix[i][i];

            }
            return 0;
        }



        public override String ToString()
        {
            string S = "";
            for (int i = 0; i < RowCount; i++)
            {
                S += "\r\n";
                for (int j = 0; j < ColumCount; j++)
                {
                    S += Matrix[i][j].ToString("F04") + "\t";
                }

                S += "\t" + Answer[i].ToString("F08");
                S += "\t" + RightPart[i].ToString("F04");
            }
            return S;
        }
    }
}
