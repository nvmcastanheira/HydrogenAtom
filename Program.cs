using MathNet.Numerics; //Documentation: https://numerics.mathdotnet.com/
using System;
using System.Collections.Generic; //https://docs.microsoft.com/en-us/dotnet/api/system.collections?view=net-5.0
using System.Linq; //https://docs.microsoft.com/en-us/dotnet/api/system.linq?view=net-5.0
using System.Numerics; //https://docs.microsoft.com/en-us/dotnet/api/system.numerics?view=net-5.0
using Microsoft.Toolkit; //https://docs.microsoft.com/en-us/windows/communitytoolkit/getting-started


//version 1.2. It runs! 
//Multiple comments are spread after each variable of importance. Remember index notation starts at 0 instead of 1 (Julia).
//You have to install some of the packages used above. Every system package is inherent. Tutorial on how to install MathNet and Microsoft.Toolkit in the doc.
//All of the important content is in the program class. There are a couple of classes/variable that I invoke in the last odd 300 lines. 
//HDF5 output is not in the code as of yet.
//Weight measurement section runs significantly faster, but time evolution runs significantly slower.
//Note: Short list of references in comments at the very end.

namespace H_atom
{
    class Program //Main program. Minimize to check other classes! References to where other classes are used can be found next to each individual class.
    {
        static void Main(string[] args)
        {
            static double laguerreCalculate(int k, double x, int alpha) //Laguerre function. Just did this recursively. Same method as code on Julia. For reference: https://en.wikipedia.org/wiki/Laguerre_polynomials, "Generalized Laguerre Polynomials" section, 4th equation down [1].

            {
                if (k == 0)
                {
                    return 1;
                }
                else if (k == 1)
                {
                    return 1 + alpha - x;
                }
                else
                {
                    return ((2 * (k - 1) + 1 + alpha - x) * laguerreCalculate((k - 1), x, alpha) - ((k - 1) + alpha) * laguerreCalculate((k - 1) - 1, x, alpha)) / k;
                }
            }

            static double legendre_Pmm(int m, double x) //Legendre polynomial of order 0. Referenced in "LegendreCalculate".
            {
                if (m == 0)
                {
                    return 1.0;
                }
                else
                {
                    double p_mm = 1.0;
                    double root_factor = Math.Sqrt(1.0 - x) * Math.Sqrt(1.0 + x);
                    double fact_coeff = 1.0;
                    int i;
                    for (i = 1; i <= m; i++)
                    {
                        p_mm *= -fact_coeff * root_factor;
                        fact_coeff += 2.0;
                    }
                    return p_mm;
                }
            }

            double legendreCalculate(int l, int m, double x) //Code for order != 0 Legendre polynomials. Made to replicate GNU special functions library: https://www.gnu.org/software/gsl/doc/html/specfunc.html#legendre-functions-and-spherical-harmonics. [2] For recursive Legendre relations.
            {

                double dif = l - m;
                double sum = l + m;
                double t_d = (dif == 0.0 ? 0.0 : 0.5 * dif * (Math.Log(dif) - 1.0));
                double t_s = (dif == 0.0 ? 0.0 : 0.5 * sum * (Math.Log(sum) - 1.0));
                double exp_check = 0.5 * Math.Log(2.0 * l + 1.0) + t_d - t_s;
                double epsilon = 2.2204460492503131e-16;
                double err_amp = 1.0 / (epsilon + Math.Abs(1.0 - Math.Abs(x)));
                double p_mm = legendre_Pmm(m, x);
                double p_mmp1 = x * (2 * m + 1) * p_mm;

                if (l == m)
                {
                    double val = p_mm;
                    double err = err_amp * 2.0 * epsilon * Math.Abs(p_mm);
                    return val;
                }
                else if (l == m + 1)
                {
                    double val = p_mmp1;
                    double err = err_amp * 2.0 * epsilon * Math.Abs(p_mmp1);
                    return val;
                }
                else
                {
                    /* upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
                     * start at P(m,m), P(m+1,m)
                     */

                    double p_ellm2 = p_mm;
                    double p_ellm1 = p_mmp1;
                    double p_ell = 0.0;
                    int ell;

                    for (ell = m + 2; ell <= l; ell++)
                    {
                        p_ell = (x * (2 * ell - 1) * p_ellm1 - (ell + m - 1) * p_ellm2) / (ell - m);
                        p_ellm2 = p_ellm1;
                        p_ellm1 = p_ell;
                    }

                    double val = p_ell;
                    double err = err_amp * (0.5 * (l - m) + 1.0) * epsilon * Math.Abs(p_ell);

                    return val;
                }

            } //Legendre function derived from the GSL code. Returns the exact values of the Legendre function. Documentation on "sf_legendre_Plm" is available on aforementioned link.            

            double get_psi_r(double r, int n, int l, int m) //Returns radial Function. [3]
            {
                double rho1 = 2.0 * r / n;
                double rho = Math.Max(rho1, Math.Pow(10, -9));
                double first = Math.Pow(2.0 / n, 3.0) / (2 * n);
                double second = SpecialFunctions.Factorial(n - l - 1) / Math.Pow(SpecialFunctions.Factorial(n + l), 3);
                double fn1 = Math.Abs(first * second);
                double my_norm = Math.Sqrt(fn1);
                double lg_fn1 = laguerreCalculate(n - l - 1, rho, 2 * l + 1) * Math.Exp(-rho / 2) * Math.Pow(rho, l) * my_norm;
                return lg_fn1; //Returns equivalent result to get_psi_r in Julia code. Consult [3] for wavefunction for the Hydrogen atom.
            }

            double get_psi_th(double cth, int l, int m) //Returns Spherical Harmonics. [3]
            {
                double first = ((((2 * l + 1) / (4 * Math.PI) * SpecialFunctions.Factorial(l - Math.Abs(m)) / SpecialFunctions.Factorial(l + Math.Abs(m)))));
                double my_norm = Math.Pow(first, 0.5);
                double lg_fn2 = legendreCalculate(l, Math.Abs(m), Math.Cos(cth));
                double result = lg_fn2 * my_norm;
                if (m >= 0)
                {
                    return Math.Pow(-1, m) * result;
                }
                else
                {
                    return result; //Returns equivalent result to get_psi_th in Julia code. Consult [3] for wavefunction for the Hydrogen atom.
                }
            }

            Complex i = new Complex(0, 1); //Complex number 0x+i. Documentation for complex numbers: https://docs.microsoft.com/en-us/dotnet/api/system.numerics.complex?view=net-5.0.
            double dt = 0.2;
            double Gamma_decay = 1.5;  //Overall scale of decay//
            //double p_meas_each = 0.1; //Probability to measure during given time step//
            double sigma_meas = 1.0; //Uncertainty in position for given measurement//
            double n_max = 3;
            double dx = 0.5; //Voxel grid step//
            double x_max = 10.25; //Voxel grid max - should not be multiple of dx//
            double r_meas = 5.0;
            double t_max = 10.0; //Evolve to max time for now//
            //double rand_seed = 0; Note: This is likely to be used for SSE once that is integrated.
            //Here, I just wrote the hydrogen parameters instead of opening up an HDF file. Variables in comments are not currently in use.

            Complex get_psi(double r, double cth, double phi, int n, int l, int m) //Full Psi composed of solutions e^i*m*phi*R(r)*O(theta,phi) to the schrodinger equation for Hydrogen [3]. 
            {
                return Complex.Exp(i * m * phi) * get_psi_r(r, n, l, m) * get_psi_th(cth, l, m);
            }

            List<double> x1 = new List<double>(); //Empty list to add x,y,z values. You'll see this a lot. Only lists have an equivalent append method, so lists will be used for some variables and may be converted to arrays later. 
            List<double> y1 = new List<double>();
            List<double> z1 = new List<double>();

            //Voxel lists for x,y,z. Maps every single permutation of x,y,z for the voxels of interest from -x_max to x_max+dx
            for (double k1 = -x_max; k1 <= x_max + dx / 2; k1 += dx) //Documentation introducing for statements in C#: https://docs.microsoft.com/en-us/dotnet/csharp/language-reference/keywords/for. 
            {
                for (double k2 = -x_max; k2 <= x_max + dx / 2; k2 += dx)
                {
                    for (double k3 = -x_max; k3 <= x_max + dx / 2; k3 += dx)
                    {
                        x1.Add(k1);
                        y1.Add(k2);
                        z1.Add(k3);
                    }
                }
            }

            //Voxel list for r.
            List<double> r1 = new List<double>();

            for (int w = 0; w < x1.Count; w++)
            {
                r1.Add(Math.Sqrt(Math.Pow(x1[w], 2) + Math.Pow(y1[w], 2) + Math.Pow(z1[w], 2)));
            }

            //Voxel list for theta.
            List<double> cth1 = new List<double>();

            for (int w = 0; w < x1.Count; w++)
            {
                cth1.Add(z1[w] / r1[w]);
            }

            List<double> phi1 = new List<double>();

            //Voxel list for phi (azimuthal angle).
            for (int w = 0; w < x1.Count; w++)
            {
                phi1.Add(Math.Atan2(y1[w], x1[w])); //Here, I use the Atan2 method. Atan1 is with regards to (1, y), and would have a phase difference when compared to the julia code.
            }

            //All the above are unsorted.

            var sorted = r1 //Lambda code; output explained below. Documentation for lambda expressions: https://docs.microsoft.com/en-us/dotnet/csharp/language-reference/operators/lambda-expressions.
                .Select((x, i) => new KeyValuePair<double, int>(x, i))
                .OrderBy(x => x.Key)
                .ToList();

            List<double> R = sorted.Select(x => x.Key).ToList(); //Here, the sorted r values are returned in a list.
            List<int> inds = sorted.Select(x => x.Value).ToList(); //Here, the INDICES of the sorted r values are returned in a list. 
            //Note: We will sort every voxel by r indices (lowest to highest). 

            //Lists that will have voxels sorted by r indices.
            List<double> X = new List<double>();
            List<double> Y = new List<double>();
            List<double> Z = new List<double>();
            List<double> CTH = new List<double>();
            List<double> PHI = new List<double>();

            foreach (int b in inds) //Sorting every variable by r indices.
            {
                X.Add(x1[b]);
                Y.Add(y1[b]);
                Z.Add(z1[b]);
                CTH.Add(cth1[b]);
                PHI.Add(phi1[b]);
            }

            int grid_sz = 0; //Grid size for r, x, y, z, theta, phi.

            for (int j = 0; j < R.Count; j++)
            {
                if (R[j] > x_max)
                {
                    grid_sz = j; //grid_sz = the index of the first r element that is higher than x_max.
                    break;
                }
            }


            //Here, we take the values [1:grid_sz] for the sorted lists, and output a list with only these values.//
            var x2 = X.Take(grid_sz); //Take method outputs every value until the argument.
            List<double> x = x2.ToList();
            var y2 = Y.Take(grid_sz);
            List<double> y = y2.ToList();
            var z2 = Z.Take(grid_sz);
            List<double> z = z2.ToList();
            var r2 = R.Take(grid_sz);
            List<double> r = r2.ToList();
            var cth2 = CTH.Take(grid_sz);
            List<double> cth = cth2.ToList();
            var phi2 = PHI.Take(grid_sz);
            List<double> phi = phi2.ToList();


            //n,l,m list, with n_max defined prior.//
            List<int> n = new List<int>();
            List<int> l = new List<int>();
            List<int> m = new List<int>();

            for (int nval = 1; nval <= n_max; nval++)
            {
                for (int lval = 0; lval < nval; lval++)
                {
                    for (int mval = -lval; mval <= lval; mval++)
                    {
                        n.Add(nval);
                        l.Add(lval);
                        m.Add(mval);
                    }
                }
            }
            //The size of each list is equal to the sum from i=1 to i=n_max of (i+1)(i+2)/2. I.e # of permutations of n,l,m.           

            Console.WriteLine("Initializing matrices...");
            Console.WriteLine("Determining energy differences...");

            List<double> E = new List<double>();
            //Energy level calculations//
            for (int w = 0; w < n.Count; w++)
            {
                E.Add(-1 * Math.Pow(1 / Convert.ToDouble(n[w]), 2)); //If we wanted the actual energy for the hydrogen atom, you would multiply these by the Rydberg constant: 1.09678*10^7 m^-1.
            }

            int basis_sz = n.Count; //E.g: basis_sz = 14 for n_max = 3//

            Complex[,] rho = new Complex[basis_sz, basis_sz]; //Empty array representing the initial value of rho. Dimensions are stated, and values are 0 unless otherwise stated.
            rho[0, 0] = 1; //Ground state matrix. Play around with multiple states. Make sure tr(rho) = 1.
            Complex[,] E_diff = new Complex[basis_sz, basis_sz]; //Array with energy differences (in Rydbergs).

            for (int row = 0; row < basis_sz; row++)
            {
                for (int col = 0; col < basis_sz; col++)
                {
                    E_diff[row, col] = E[row] - E[col];
                }
            }

            Console.WriteLine("Reading in Decay Matrices..."); //Gamma values start to be used here. Refer to the "variables" class if needed.


            //Here, I simply imported the hdf5 matrices by hand, because of issues with hdf5 methods. I.e, I copy-pasted the matrices that you get from the file. Works just fine, but feel free to import them with hdf5 if desired.//
            List<double> Gamma_lst = variables.Gamma_lst0.ToList(); //Values of Gamma attained from other code.
            List<int> gamma_ind_lst = variables.gamma_ind_lst0.ToList(); //Indices of the Gamma values.
            List<int> delta_ind_lst = variables.delta_ind_lst0.ToList(); //Indices of the delta values.

            double max_rate = 0; //Create max_rate double.

            for (int j = 0; j < gamma_ind_lst.Count; j++)
            {
                if (gamma_ind_lst[j] == 1)
                {
                    max_rate = Math.Max(max_rate, Gamma_lst[j]); //Returns max rate, which corresponds to the first index of the Gamma list.
                }

            }

            Gamma_lst.Add(0.5 * max_rate); //Adds half the max rate to the Gamma list.
            gamma_ind_lst.Add(1); //Adds index 1 to the Gamma indices.
            delta_ind_lst.Add(2); //Adds index 2 to the Delta indices.

            string s = (max_rate * 0.5) + "";

            Console.WriteLine("Adding artificial decay from 2,0,0 -> 1,0,0 at rate " + s);

            double[,] decay_matrix_coh = new double[basis_sz, basis_sz]; //Decay matrix for off-diagonal terms. I.e, Gamma values.
            double[,] decay_matrix_ampl = new double[basis_sz, basis_sz]; //Decay matrix for diagonal terms.

            for (int j = 0; j < gamma_ind_lst.Count; j++) //Determines the values of the decay matrices for time evolution. 
            {
                int gi = gamma_ind_lst[j] - 1;
                int di = delta_ind_lst[j] - 1;
                double G = Gamma_lst[j];
                if (gi >= basis_sz || di >= basis_sz) //Allows only index values gi < basis_sz and di < basis_sz.
                {
                    continue;
                }
                decay_matrix_ampl[gi, di] += G; //Adds gamma values to non-diagonal entries.
                decay_matrix_ampl[di, di] -= G; //Subtracts gamma values for diagonal entries.
                for (int alpha = 0; alpha < basis_sz; alpha++)
                {
                    if (alpha != di)
                    {
                        decay_matrix_coh[di, alpha] -= G / 2;
                        decay_matrix_coh[alpha, di] -= G / 2;
                    } //Subtracts gamma/2 for non-diagonal entries of this matrix. Diagonal values of this matrix are always 0.
                }
            }

            double[,] exp_decay_matrix_ampl = new double[basis_sz, basis_sz]; //New list. Used during time evolution.

            for (int j = 0; j < basis_sz; j++)
            {
                for (int k = 0; k < basis_sz; k++)
                {
                    exp_decay_matrix_ampl[j, k] = Math.Exp(dt * decay_matrix_ampl[j, k] * Gamma_decay); //Simply the decay matrix in an exp function. Will be used for diagonal terms in time evolution.
                }
            }

            Console.WriteLine("Constructing wfns at each site...");

            Complex[][] psi_grid0 = new Complex[basis_sz][]; //This will be our first psi_grid. The final psi_grid is the array below. It is a jagged array instead of a 2-d array, simply for convenience. Documentation on jagged arrays: https://docs.microsoft.com/en-us/dotnet/csharp/programming-guide/arrays/jagged-arrays.
            Complex[][] psi_grid = new Complex[basis_sz][];//grid_sz is the second parameter. This is equivalent to the first psi_grid seen on the Julia code. 

            for (int j = 0; j < basis_sz; j++) //This code defines our second dimension for the jagged arrays: so that we have a basis_sz X grid_sz dimensional array.
            {
                psi_grid0[j] = new Complex[grid_sz];
                psi_grid[j] = new Complex[grid_sz];
            }

            for (int j = 0; j < basis_sz; j++)
            {
                if (j % 10 == 1)
                {
                    Console.WriteLine("At j = " + j + " out of " + basis_sz);
                }
                for (int k = 0; k < grid_sz; k++)
                {
                    psi_grid0[j][k] = get_psi(r[k], cth[k], phi[k], n[j], l[j], m[j]); //We assign psi values to this grid. Moving horizontally (column-wise) on the array represents a change in r, cth, phi. Moving vertically represents a change in quantum numbers.
                }
            }

            for (int j = 0; j < grid_sz; j++) //Code to normalize the psi grid.
            {
                if (j % 10 == 0)
                {
                    Console.WriteLine("At j = " + j + " out of " + grid_sz);
                }
                Complex[] psicol1 = new Complex[basis_sz]; //Here I need to invoke a couple more lists to make the process a bit easier to work with.
                Complex[] psicol3 = new Complex[basis_sz];

                var psicolr = ArrayExtensions.GetColumn<Complex>(psi_grid0, j); //Part of Microsoft.Tools. Gets jth column of psi_grid0, and maps to var "psicolr".
                Complex[] psicol = psicolr.ToArray(); //Turns psicolr to array.

                double[] psicol2 = new double[basis_sz];

                for (int b = 0; b < basis_sz; b++) //Operations from here to the end of the loop are such that: psi_grid[:,j]=psi_grid[:,j]./(sum((abs.(psi_grid[:,j])).^2))^0.5
                {
                    psicol1[b] = Complex.Abs(psicol[b]);
                    psicol1[b] = Complex.Pow(psicol1[b], 2);
                    psicol2[b] = psicol1[b].Real;
                }
                double sum = psicol2.Sum();

                double Normal = Math.Pow(sum, 0.5);

                for (int b = 0; b < basis_sz; b++)
                {
                    psicol3[b] = psicol[b] / Normal; //Normalizing each column.
                    psi_grid[b][j] = psicol3[b]; //Map every value of jth column to row b, column j of psi_grid.
                }
            }

            Console.WriteLine("Constructing meas weights..."); //This part runs rather fast compared to the Julia code (no idea why). weight_meas isn't currently used in the code.

            double[,] weight_meas = new double[grid_sz, grid_sz]; //Weight measurements. Not used as of yet.
            for (int row = 0; row < grid_sz; row++) //Maps weight measurement values.
            {
                if (row % 10 == 0)
                {
                    Console.WriteLine("At row = " + row + " out of " + grid_sz);
                }
                for (int col = 0; col < grid_sz; col++)
                {
                    weight_meas[row, col] = Math.Exp(-(Math.Pow((x[row] - x[col]), 2) + Math.Pow((y[row] - y[col]), 2) + Math.Pow((z[row] - z[col]), 2) / (2 * Math.Pow(sigma_meas, 2))));
                }
            }

            Console.WriteLine("Creating measurement operators from weights + wfns...");

            Complex[,,] M_meas = new Complex[1, basis_sz, basis_sz]; //Essentially a basis_sz X basis_sz matrix. First dimension is 1 for use later on.

            int g0 = 0;

            for (int g = 0; g < grid_sz; g++)
            {
                if (r[g] >= r_meas) //Finds when r is greater than r_meas, and maps its index to g0.
                {
                    g0 = g;
                    break;
                }

            }

            for (int g = 0; g < grid_sz; g++) //Assigns measurement operators for time evolution.
            {
                if (g % 10 == 0)
                {
                    Console.WriteLine("At g = " + g + " out of " + grid_sz);
                }

                for (int j = 0; j < basis_sz; j++)
                {
                    for (int k = 0; k < basis_sz; k++)
                    {

                        M_meas[0, j, k] += Math.Exp(-(Math.Pow(x[g] - x[g0], 2) + Math.Pow(y[g] - y[g0], 2) + Math.Pow(z[g] - z[g0], 2)) / 4.0) * psi_grid[j][g] * Complex.Conjugate(psi_grid[k][g]);

                    }
                }
            }

            List<double> t = new List<double>(); //Time list.

            for (double j = 0; j <= t_max + (dt / 2); j += dt) //Time values, in intervals of dt.
            {
                t.Add(j);
            }

            double[,] p = new double[grid_sz, t.Count]; //Probability for each voxel. From the dimensions, its clear that we are measuring the probability per voxel for the given time value.
            double[,] p_basis = new double[basis_sz, t.Count]; //Each column of this array is the diagonal of the rho matrix for a given time t.
            double[] meas_ind = new double[t.Count]; //Never referenced.
            double[] p_sum = new double[grid_sz]; //Used in initial observable code.

            Console.WriteLine("Evaluating observable for initial state...");

            for (int g = 0; g < grid_sz; g++) //Code for probability prior to time evolution (t=0).
            {
                var psicolw = ArrayExtensions.GetColumn<Complex>(psi_grid, g); //Maps each column of psi_grid to a vector.
                Complex[,] psicol4 = new Complex[basis_sz, 1]; //Placeholder array for operations. Cannot use method "TransposeRowsAndColumns" if we use psicol5, because the method only works with two-dimensional arrays.
                Complex[] psicol5 = psicolw.ToArray(); //Turns var psicolw to a 1D array.

                for (int j = 0; j < basis_sz; j++)
                {
                    psicol4[j, 0] = psicol5[j];
                    p_basis[j, 0] = rho[j, j].Real;
                }
                Complex[,] first = matrix_mult_C.Multiply(rho, psicol4);

                Complex[,] second = matrix_mult_C.Multiply(TransposeRowsColumnsExtension.TransposeRowsAndColumns(psicol4), first);

                p[g, 0] = second[0, 0].Real; //This and prior operations are equivalent to: p[g,1]=real(psi_grid[:,g]'*(rho*psi_grid[:,g])).
                p_sum[g] = p[g, 0]; //Maps all values of the probability voxels to a vector for normalization.
            }

            for (int g = 0; g < grid_sz; g++) //Normalizes probability voxels for time t = 0.
            {
                p[g, 0] = p[g, 0] / p_sum.Sum();
            }

            Console.WriteLine("Starting time evolution");

            //Variables are instantiated before time evolution code so they are not overwritten during multiple iterations.
            //Most arrays are placeholders for operations. 
            //Will go into more detail for this section, since this is where the code differs the most between the Julia code (way more loops/variables are necessary for slicing and other operations).

            Complex[,] rho_inds = new Complex[basis_sz, 1]; //Rho diagonal values.
            Complex[,] comp_exp = new Complex[basis_sz, basis_sz]; //Maps exp_decay_matrix_ampl to a complex value matrix. This is necessary so that variable types don't mix.
            Complex[,] decay_rho = new Complex[basis_sz, basis_sz]; //Matrix multiplication array of comp_exp * rho_inds.
            Complex[,] M_meas1 = new Complex[basis_sz, basis_sz]; //Same as M_meas, but on a 2D array instead.
            Complex[,] psicol_f = new Complex[basis_sz, 1]; //Array for each column of psi_grid. Needed for matrix_mult.
            Complex[] tr_rho = new Complex[basis_sz]; //Another array of diagonal rho values. Iterated over each time t.
            double[] p_sum1 = new double[grid_sz]; //List used for normalizing p.
            double[] p_after = new double[grid_sz]; //Probability voxels after time evolution. I.e, voxels at the last t-value.            
            Complex trace = new Complex(0, 0); //Trace of rho.

            for (int row = 0; row < basis_sz; row++)
            {
                for (int col = 0; col < basis_sz; col++)
                {
                    comp_exp[row, col] = exp_decay_matrix_ampl[row, col]; //Map each exp_decay_matrix_ampl element to 2D array comp_exp.
                }
            }

            for (int j = 1; j < t.Count; j++) //Time evolution from t = dt to t = t_max
            {
                if (j % 5 == 0)
                {
                    Console.WriteLine("At step j = " + j + " out of " + t.Count);
                }

                for (int row = 0; row < basis_sz; row++)
                {
                    for (int col = 0; col < basis_sz; col++)
                    {
                        rho[row, col] = Complex.Exp((Complex.Add(-i * E_diff[row, col], Gamma_decay * decay_matrix_coh[row, col])) * dt) * rho[row, col]; //Solving for the time-evolved density matrix.
                    }
                }

                for (int row = 0; row < basis_sz; row++)
                {
                    rho_inds[row, 0] = rho[row, row]; //Map diagonal terms of rho to array rho_inds.                   
                }

                decay_rho = matrix_mult_C.Multiply(comp_exp, rho_inds); //Multiplies comp_exp, rho_inds (both are complex matrices). 

                for (int row = 0; row < basis_sz; row++)
                {
                    rho[row, row] = decay_rho[row, 0]; //decay_rho elements are mapped to diagonal of rho. When combined with lines 543-547, this is the same as: rho[dinds]=exp_decay_matrix_ampl*rho[dinds].
                }

                for (int g = 0; g < grid_sz; g++)
                {
                    var psicolw = ArrayExtensions.GetColumn<Complex>(psi_grid, g);
                    Complex[] psicol5 = psicolw.ToArray(); //552-553 make an array for each column of psi_grid.

                    for (int col = 0; col < basis_sz; col++)
                    {
                        psicol_f[col, 0] = psicol5[col];
                    }

                    Complex[,] third = matrix_mult_C.Multiply(rho, psicol_f);
                    Complex[,] fourth = matrix_mult_C.Multiply(psicol_f.TransposeRowsAndColumns(), third);

                    p[g, j] = fourth[0, 0].Real; //560-563 is the same as: p[g,j]=real(psi_grid[:,g]'*(rho*psi_grid[:,g]))

                    p_sum1[g] = p[g, j]; //Maps every element of each jth column of p to normalize the matrix.
                }

                for (int g = 0; g < grid_sz; g++)
                {
                    p[g, j] /= p_sum1.Sum(); //Normalizes each column of the p.
                }

                if (j == 1) //Measure at step 1.
                {
                    for (int row = 0; row < basis_sz; row++)
                    {
                        for (int col = 0; col < basis_sz; col++)
                        {
                            M_meas1[row, col] = M_meas[0, row, col]; //Maps each value of M_meas to a 2D array M_meas1. Important for matrix multiplication.
                        }
                    }

                    Complex[,] fifth = matrix_mult_C.Multiply(rho, M_meas1);
                    Complex[,] sixth = matrix_mult_C.Multiply(M_meas1, fifth);

                    for (int row = 0; row < basis_sz; row++)
                    {
                        for (int col = 0; col < basis_sz; col++)
                        {
                            rho[row, col] = sixth[row, col]; //597-606 is the same as: rho=M_meas[1,:,:]*rho*M_meas[1,:,:]
                        }
                    }

                    for (int row = 0; row < basis_sz; row++)
                    {
                        trace += rho[row, row]; //Trace of rho.
                    }

                    for (int row = 0; row < basis_sz; row++)
                    {
                        for (int col = 0; col < basis_sz; col++)
                        {
                            rho[row, col] /= trace; //Normalize rho: rho = rho/tr(rho)
                        }
                    }
                }

                for (int g = 0; g < grid_sz; g++)
                {
                    var psicolw = ArrayExtensions.GetColumn<Complex>(psi_grid, g);
                    Complex[] psicol5 = psicolw.ToArray();

                    for (int row = 0; row < basis_sz; row++)
                    {
                        psicol_f[row, 0] = psicol5[row];
                    }
                    Complex[,] third = matrix_mult_C.Multiply(rho, psicol_f);
                    Complex[,] fourth = matrix_mult_C.Multiply(psicol_f.TransposeRowsAndColumns(), third);

                    p_after[g] = fourth[0, 0].Real;

                }

                for (int g = 0; g < grid_sz; g++)
                {
                    p_after[g] = p_after[g] / p_after.Sum();
                }
                for (int col = 0; col < basis_sz; col++)
                {
                    p_basis[col, j] = rho[col, col].Real;
                }

            }


            for (int w = 0; w < 10; w++)
            {
                Console.WriteLine("p_after[" + w + "] = " + p_after[w]);
            }

            for (int w = 0; w < 4; w++)
            {
                for (int ip = 0; ip < 4; ip++)
                {
                    Console.WriteLine("p[" + w + "][" + ip + "] = " + p[w, ip]);
                    Console.WriteLine("pbasis[" + w + "][" + ip + "] = " + p_basis[w, ip]);

                }
            }

            /*# Then save data
            fout = h5open("data_hydrogen.h5", "w")
            write(fout, "p", p)
            write(fout, "p_basis", p_basis)
            write(fout, "meas_ind", meas_ind)
            write(fout, "x", convert(Array{ Float64,1},x))
            write(fout, "y", convert(Array{ Float64,1},y))
            write(fout, "z", convert(Array{ Float64,1},z))
            write(fout, "n", convert(Array{ Int64,1},n))
            write(fout, "l", convert(Array{ Int64,1},l))
            write(fout, "m", convert(Array{ Int64,1},m))
            write(fout, "t", collect(t))
            close(fout)
            This is how the hdf5 output code looks in Julia. You will likely need to write some code to output a similar hdf5 file. */

            //and that ends run_main.
        }
    }
    public static class variables
    {
        public static double[] Gamma_lst0 = {0.27746455641634377, 0.5549291128326874, 0.27746455641634377, 0.04449479088272272, 0.08898958176544547, 0.04449479088272272, 0.015461882472182085, 0.030923764944364173, 0.015461882472182085, 0.007259554259128699, 0.014519108518257402, 0.007259554259128699, 0.0040111677909952246, 0.008022335581990449, 0.0040111677909952246, 0.0024571208009544114, 0.004914241601908823, 0.0024571208009544114, 0.0016169789990632545, 0.003233957998126509, 0.0016169789990632545, 0.0011219051096993315, 0.0022438102193986635, 0.0011219051096993315, 0.0008107916025319948, 0.0016215832050639899, 0.0008107916025319948, 4.4999998473627345, 8.99999969472547, 4.4999998473627345, 1.5655155833651695, 3.1310311667303394, 1.5655155833651695, 0.2740390078918631, 0.5480780157837262, 0.2740390078918631, 0.09984666776405761, 0.19969333552811525, 0.09984666776405761, 0.04866602261037278, 0.09733204522074554, 0.04866602261037278, 0.02774572255578744, 0.05549144511157489, 0.02774572255578744, 0.017451862377448536, 0.034903724754897066, 0.017451862377448536, 0.011746225317491178, 0.023492450634982356, 0.011746225317491178, 0.008308856731346707, 0.016617713462693414, 0.008308856731346707, 0.1467670865488975, 4.508688681703181, 4.508688681703185, 0.7514481136171964, 0.02435902506290565, 0.5846248786214672, 0.5846248786214676, 0.09743747977024449, 0.008666103020546342, 0.19015907500193455, 0.19015907500193469, 0.03169317916698908, 0.004171371459318171, 0.08759880838851812, 0.08759880838851819, 0.014599801398086348, 0.0023603132031175767, 0.048339214457524426, 0.04833921445752447, 0.008056535742920735, 0.0014774058571157148, 0.02978449795019917, 0.0297844979501992,
0.00496408299169986, 0.00099108761933465, 0.019770261626456172, 0.01977026162645619, 0.0032950436044093603, 0.0006993987475717612, 0.01384808930377602, 0.01384808930377603, 0.002308014883962669, 0.2935341730977951, 2.254344340851592, 6.011584908937574, 2.254344340851592, 0.048718050125811305, 0.2923124393107337, 0.7794998381619562, 0.2923124393107337, 0.017332206041092684, 0.09507953750096731, 0.2535454333359127, 0.09507953750096731, 0.008342742918636344, 0.04379940419425909, 0.11679841118469081, 0.04379940419425909, 0.0047206264062351535, 0.024169607228762227, 0.06445228594336591, 0.024169607228762227, 0.00295481171423143, 0.014892248975099593, 0.039712663933598895, 0.014892248975099593, 0.0019821752386693, 0.009885130813228091, 0.02636034883527489, 0.009885130813228091, 0.0013987974951435223, 0.006924044651888015, 0.01846411907170136, 0.006924044651888015, 0.1467670865488975, 0.7514481136171964, 4.508688681703185, 4.508688681703181, 0.02435902506290565, 0.09743747977024449, 0.5846248786214676, 0.5846248786214672, 0.008666103020546342, 0.03169317916698908, 0.19015907500193469, 0.19015907500193455, 0.004171371459318171, 0.014599801398086348, 0.08759880838851819, 0.08759880838851812, 0.0023603132031175767, 0.008056535742920735, 0.04833921445752447, 0.048339214457524426, 0.0014774058571157148, 0.00496408299169986, 0.0297844979501992, 0.02978449795019917, 0.00099108761933465, 0.0032950436044093603, 0.01977026162645619, 0.019770261626456172, 0.0006993987475717612, 0.002308014883962669, 0.01384808930377603, 0.01384808930377602, 27.000000614395592, 54.00000122879118, 27.000000614395592, 4.985605465757618, 9.971210931515238, 4.985605465757618, 0.8509468502371562, 1.7018937004743124, 0.8509468502371562, 0.3083644038923534, 0.6167288077847067, 0.3083644038923534, 0.15086877815703098, 0.301737556314062, 0.15086877815703098, 0.08665789681053229, 0.1733157936210646, 0.08665789681053229, 0.05499372665351863, 0.10998745330703726, 0.05499372665351863, 0.03735647448688872, 0.07471294897377743, 0.03735647448688872, 20.25000003016822, 20.250000030168238, 3.3750000050280353, 0.995142952684912, 11.447088155234143, 11.447088155234152, 1.9078480258723565, 0.15669066648125038, 1.7621872133064354, 1.7621872133064371, 0.29369786888440574, 0.05458225175756314, 0.6062690305513442, 0.6062690305513446, 0.10104483842522398, 0.02610871252664321, 0.2878375238576533, 0.28783752385765354, 0.047972920642942195, 0.014785384696582601, 0.16221412718328634, 0.16221412718328648, 0.027035687863881053, 0.009293929673655865, 0.10162833011764827, 0.10162833011764838, 0.016938055019608038, 0.006270953336060207, 0.06841012128265872, 0.06841012128265879, 0.011401686880443119,
10.125000015084117, 27.000000040224297, 10.125000015084117, 1.9902859053698239, 5.723544077617074, 15.262784206978855, 5.723544077617074, 0.31338133296250087, 0.8810936066532183, 2.3495829510752473, 0.8810936066532183, 0.1091645035151263, 0.3031345152756722, 0.8083587074017922, 0.3031345152756722, 0.052217425053286415, 0.14391876192882672, 0.3837833651435377, 0.14391876192882672, 0.029570769393165202, 0.08110706359164323, 0.21628550291104848, 0.08110706359164323, 0.01858785934731173, 0.05081416505882417, 0.13550444015686436, 0.05081416505882417, 0.012541906672120415, 0.034205060641329395, 0.09121349504354498, 0.034205060641329395, 3.3750000050280353, 20.250000030168238, 20.25000003016822, 0.995142952684912, 1.9078480258723565, 11.447088155234152, 11.447088155234143, 0.15669066648125038, 0.29369786888440574, 1.7621872133064371, 1.7621872133064354, 0.05458225175756314, 0.10104483842522398, 0.6062690305513446, 0.6062690305513442, 0.02610871252664321, 0.047972920642942195, 0.28783752385765354, 0.2878375238576533, 0.014785384696582601, 0.027035687863881053, 0.16221412718328648, 0.16221412718328634, 0.009293929673655865,
0.016938055019608038, 0.10162833011764838, 0.10162833011764827, 0.006270953336060207, 0.011401686880443119, 0.06841012128265879, 0.06841012128265872, 0.3391729826482571, 22.426948052586322, 14.951298701724218, 1.4951298701724196, 0.046618710619140864, 2.360072144574329, 1.57338142971622, 0.15733814297162177, 0.01515672445499474, 0.6929042174502217, 0.4619361449668146, 0.046193614496681384, 0.006964955598796279, 0.3022293442173216, 0.2014862294782145, 0.02014862294782142, 0.003845074190604645, 0.16179274158623497, 0.1078618277241567, 0.010786182772415652, 0.0023757307014370967, 0.09799871076539773, 0.06533247384359851, 0.0065332473843598415, 0.0015835698727462483, 0.06443651206688397, 0.04295767471125599, 0.004295767471125592, 0.3391729826482575, 0.1695864913241287, 14.951298701724218, 23.922077922758774, 4.485389610517265, 0.046618710619140906, 0.023309355309570446, 1.57338142971622, 2.5174102875459545, 0.47201442891486595, 0.015156724454994754, 0.007578362227497375, 0.4619361449668146, 0.739097831946904, 0.13858084349004435, 0.0069649555987962835, 0.0034824777993981413, 0.2014862294782145, 0.3223779671651436, 0.06044586884346434, 0.0038450741906046475, 0.0019225370953023237, 0.1078618277241567, 0.1725789243586509, 0.03235854831724701, 0.002375730701437099, 0.0011878653507185492, 0.06533247384359851, 0.10453195814975771, 0.019599742153079545, 0.0015835698727462494, 0.0007917849363731246, 0.04295767471125599, 0.06873227953800964, 0.012887302413376793, 0.0565288304413762, 0.45223064353100956, 0.0565288304413762, 8.970779221034533, 26.91233766310356, 8.970779221034537, 0.007769785103190146, 0.06215828082552116, 0.007769785103190146, 0.9440288578297323, 2.8320865734891933, 0.9440288578297328, 0.0025261207424991236, 0.02020896593999298, 0.0025261207424991236, 0.2771616869800888, 0.8314850609402654, 0.277161686980089, 0.0011608259331327133, 0.009286607465061705, 0.0011608259331327133, 0.12089173768692878, 0.3626752130607858, 0.1208917376869288, 0.0006408456984341074, 0.0051267655874728596, 0.0006408456984341074, 0.06471709663449406, 0.19415128990348185, 0.06471709663449406, 0.0003959551169061829, 0.0031676409352494627, 0.0003959551169061829, 0.03919948430615912, 0.11759845291847719, 0.039199484306159126, 0.00026392831212437477, 0.0021114264969949977, 0.00026392831212437477, 0.025774604826753603, 0.07732381448026071, 0.025774604826753613, 0.1695864913241287, 0.3391729826482575, 4.485389610517267, 23.922077922758774, 14.951298701724221, 0.023309355309570446, 0.046618710619140906, 0.4720144289148661, 2.5174102875459545, 1.5733814297162205, 0.007578362227497375, 0.015156724454994754, 0.13858084349004438, 0.739097831946904, 0.46193614496681473, 0.0034824777993981413, 0.0069649555987962835, 0.06044586884346436, 0.3223779671651436, 0.20148622947821457, 0.0019225370953023237, 0.0038450741906046475, 0.032358548317247014, 0.1725789243586509, 0.10786182772415674, 0.0011878653507185492, 0.002375730701437099, 0.019599742153079552, 0.10453195814975771,
0.06533247384359851, 0.0007917849363731246, 0.0015835698727462494, 0.0128873024133768, 0.06873227953800964, 0.042957674711256, 0.3391729826482571, 1.4951298701724196, 14.951298701724218, 22.42694805258631,
0.046618710619140864, 0.15733814297162177, 1.57338142971622, 2.3600721445743287, 0.01515672445499474, 0.046193614496681384, 0.4619361449668146, 0.6929042174502215, 0.006964955598796279, 0.02014862294782142, 0.2014862294782145, 0.3022293442173216, 0.003845074190604645, 0.010786182772415652, 0.1078618277241567, 0.16179274158623494, 0.0023757307014370967, 0.0065332473843598415, 0.06533247384359851, 0.0979987107653977, 0.0015835698727462483, 0.004295767471125592, 0.04295767471125599, 0.06443651206688394, 90.00000050774808, 180.00000101549622, 90.00000050774808, 12.092228238868211, 24.184456477736415, 12.092228238868211, 1.9889303840051855, 3.9778607680103715, 1.9889303840051855, 0.7085440008570948, 1.4170880017141896, 0.7085440008570948, 0.34421233982135824, 0.6884246796427166, 0.34421233982135824, 0.1973622995843607, 0.3947245991687214, 0.1973622995843607, 0.12538853790382562, 0.25077707580765124, 0.12538853790382562, 86.40000012651983, 86.4000001265199, 14.400000021086635, 3.527092307755915, 24.37165365138013, 24.371653651380146, 4.061942275230019, 0.5329742115195469, 3.8481990380511237, 3.8481990380511264, 0.6413665063418537, 0.18154240204297678, 1.3373642088129594, 1.3373642088129605, 0.2228940348021598, 0.08583760981594031, 0.6393684290698677, 0.639368429069868, 0.10656140484497789, 0.0483509802589876, 0.3626080586518156, 0.3626080586518159, 0.060434676441969246, 0.030342140105921102, 0.22858311556113317, 0.22858311556113337, 0.038097185926855515, 43.200000063259935, 115.20000016869311, 43.200000063259935, 7.054184615511829, 12.185826825690071, 32.49553820184017, 12.185826825690071, 1.0659484230390939, 1.924099519025563, 5.13093205073483, 1.924099519025563, 0.3630848040859536, 0.6686821044064802, 1.783152278417279, 0.6686821044064802, 0.17167521963188068, 0.319684214534934, 0.8524912387598235, 0.319684214534934, 0.0967019605179752, 0.18130402932590794, 0.48347741153575413, 0.18130402932590794, 0.06068428021184221, 0.11429155778056667, 0.3047774874148443, 0.11429155778056667, 14.400000021086635, 86.4000001265199, 86.40000012651983, 3.527092307755915, 4.061942275230019, 24.371653651380146, 24.37165365138013, 0.5329742115195469, 0.6413665063418537, 3.8481990380511264, 3.8481990380511237, 0.18154240204297678, 0.2228940348021598,
1.3373642088129605, 1.3373642088129594, 0.08583760981594031, 0.10656140484497789, 0.639368429069868, 0.6393684290698677, 0.0483509802589876, 0.060434676441969246, 0.3626080586518159, 0.3626080586518156, 0.030342140105921102, 0.038097185926855515, 0.22858311556113337, 0.22858311556113317, 53.99999996793794, 35.999999978625304, 3.599999997862525, 1.8547949424693413, 42.391460566997935, 28.260973711331957, 2.8260973711331916, 0.25217932931902004, 5.744153975033745, 3.829435983355831, 0.38294359833558256, 0.08103668071385063, 1.8389642982425312, 1.2259761988283546, 0.12259761988283527, 0.03696198493077183, 0.836396939189295, 0.5575979594595302, 0.055759795945952935, 0.020328053312883676, 0.45903681351652625, 0.3060245423443509, 0.030602454234435045, 0.012544665973740234, 0.2828393205135535, 0.18855954700903574, 0.018855954700903545, 35.999999978625304, 57.59999996580055, 10.79999999358759, 1.8547949424693424, 0.9273974712346711, 28.260973711331957, 45.217557938131186, 8.478292113399586, 0.25217932931902015, 0.12608966465951008, 3.829435983355831, 6.127097573369336, 1.148830795006749, 0.0810366807138507, 0.040518340356925336, 1.2259761988283546, 1.9615619181253694, 0.36779285964850633, 0.03696198493077185, 0.018480992465385926, 0.5575979594595302, 0.8921567351352493, 0.16727938783785903, 0.02032805331288369, 0.010164026656441845, 0.3060245423443509, 0.489639267750962, 0.09180736270330524, 0.012544665973740243, 0.0062723329868701215, 0.18855954700903574, 0.30169527521445744, 0.0565678641027107, 21.599999987175188, 64.7999999615255, 21.5999999871752, 0.3091324904115569, 2.473059923292454, 0.3091324904115569, 16.95658422679918,
50.86975268039748, 16.956584226799187, 0.04202988821983667, 0.3362391057586933, 0.04202988821983667, 2.2976615900134996, 6.89298477004049, 2.2976615900135005, 0.013506113452308437, 0.10804890761846751, 0.013506113452308437, 0.735585719297013, 2.2067571578910363, 0.7355857192970132, 0.006160330821795305, 0.04928264657436244, 0.006160330821795305, 0.33455877567571823, 1.0036763270271531, 0.33455877567571835, 0.003388008885480613, 0.027104071083844904, 0.003388008885480613, 0.18361472540661064, 0.5508441762198311, 0.1836147254066107, 0.002090777662290039, 0.016726221298320312, 0.002090777662290039, 0.11313572820542146, 0.33940718461626396, 0.1131357282054215, 10.799999993587592, 57.59999996580055, 35.99999997862531, 0.9273974712346711, 1.8547949424693424, 8.478292113399588, 45.217557938131186, 28.260973711331967, 0.12608966465951008, 0.25217932931902015, 1.1488307950067496, 6.127097573369336, 3.829435983355832, 0.040518340356925336, 0.0810366807138507, 0.3677928596485065, 1.9615619181253694, 1.225976198828355, 0.018480992465385926, 0.03696198493077185, 0.1672793878378591, 0.8921567351352493, 0.5575979594595303, 0.010164026656441845, 0.02032805331288369, 0.09180736270330528, 0.489639267750962, 0.306024542344351, 0.0062723329868701215, 0.012544665973740243, 0.056567864102710724, 0.30169527521445744, 0.18855954700903577, 3.599999997862525, 35.999999978625304, 53.999999967937924, 1.8547949424693413, 2.8260973711331916, 28.260973711331957, 42.39146056699791, 0.25217932931902004, 0.38294359833558256, 3.829435983355831, 5.744153975033743, 0.08103668071385063, 0.12259761988283527, 1.2259761988283546, 1.8389642982425305, 0.03696198493077183, 0.055759795945952935, 0.5575979594595302, 0.8363969391892948, 0.020328053312883676, 0.030602454234435045, 0.3060245423443509, 0.4590368135165261, 0.012544665973740234, 0.018855954700903545, 0.18855954700903574, 0.28283932051355337, 0.5914157987683235, 69.7821973521976, 34.891098676098885, 2.492221334007058, 0.0683827854235496, 6.1271012268070235, 3.063550613403518, 0.21882504381453655, 0.02011367196214034, 1.6181234684826746, 0.8090617342413393, 0.057790123874381254, 0.008689838125755776, 0.6608568434200005, 0.3304284217100011, 0.023602030122142882, 0.004610860536255203, 0.33898103536485763, 0.16949051768242923, 0.01210646554874492, 0.002775110777797869, 0.19953391571096066, 0.09976695785548054, 0.007126211275391454, 0.3942771991788825, 0.3942771991788825, 52.336648014148345, 59.813312016169476, 7.476664002021176, 0.04558852361569973, 0.04558852361569973, 4.595325920105279, 5.251801051548885, 0.6564751314436099, 0.013409114641426897, 0.013409114641426897, 1.2135926013620093, 1.386962972985152, 0.1733703716231438, 0.0057932254171705186, 0.0057932254171705186, 0.4956426325650017, 0.56644872293143, 0.07080609036642868, 0.0030739070241701362, 0.0030739070241701362, 0.2542357765236439, 0.2905551731698784, 0.036319396646234764, 0.0018500738518652468, 0.0018500738518652468, 0.14965043678322087, 0.17102907060939507, 0.021378633826174363, 0.039427719917888196, 0.6308435186862126, 0.2365663195073295, 37.3833200101059,
74.76664002021187, 14.953328004042351, 0.004558852361569969, 0.07294163778511965, 0.02735311416941985,
3.28237565721805, 6.564751314436106, 1.3129502628872196, 0.0013409114641426882, 0.02145458342628306, 0.008045468784856142, 0.8668518581157194, 1.7337037162314404, 0.3467407432462876, 0.0005793225417170513,
0.00926916066747284, 0.0034759352503023126, 0.35403045183214343, 0.7080609036642875, 0.14161218073285733, 0.0003073907024170133, 0.004918251238672224, 0.0018443442145020825, 0.18159698323117388, 0.3631939664623481, 0.07263879329246953, 0.00018500738518652444, 0.0029601181629843975, 0.0011100443111191484, 0.10689316913087185, 0.2137863382617439, 0.042757267652348725, 0.1182831597536647, 0.7096989585219877, 0.11828315975366474, 24.92221334007059, 79.75108268822575, 24.92221334007061, 0.013676557084709919, 0.08205934250825946, 0.013676557084709925, 2.1882504381453667, 7.002401402065159, 2.188250438145368, 0.004022734392428068, 0.024136406354568388, 0.00402273439242807, 0.5779012387438129, 1.8492839639801972, 0.5779012387438133, 0.0017379676251511552, 0.010427805750906922, 0.0017379676251511559, 0.23602030122142895,
0.7552649639085709, 0.2360203012214291, 0.0009221721072510405, 0.005533032643506239, 0.0009221721072510409, 0.12106465548744924, 0.3874068975598368, 0.12106465548744931, 0.0005550221555595739, 0.003330132933357441, 0.000555022155559574, 0.07126211275391457, 0.22803876081252616, 0.0712621127539146, 0.2365663195073295, 0.6308435186862126, 0.03942771991788819, 14.953328004042355, 74.7666400202119, 37.3833200101059, 0.02735311416941985, 0.07294163778511965, 0.004558852361569967, 1.31295026288722, 6.56475131443611, 3.28237565721805, 0.008045468784856142, 0.02145458342628306, 0.0013409114641426876, 0.34674074324628773, 1.7337037162314415, 0.8668518581157194, 0.0034759352503023126, 0.00926916066747284, 0.0005793225417170509, 0.14161218073285736, 0.7080609036642881, 0.35403045183214343, 0.0018443442145020825, 0.004918251238672224, 0.0003073907024170131, 0.07263879329246956, 0.3631939664623483, 0.18159698323117388, 0.0011100443111191484, 0.0029601181629843975, 0.00018500738518652439, 0.04275726765234873, 0.21378633826174406, 0.10689316913087185, 0.39427719917888265, 0.3942771991788825, 7.476664002021176, 59.813312016169476,
52.336648014148345, 0.04558852361569976, 0.04558852361569973, 0.6564751314436099, 5.251801051548885, 4.595325920105279, 0.013409114641426904, 0.013409114641426897, 0.1733703716231438, 1.386962972985152, 1.2135926013620093, 0.005793225417170521, 0.0057932254171705186, 0.07080609036642868, 0.56644872293143, 0.4956426325650017, 0.003073907024170138, 0.0030739070241701362, 0.036319396646234764, 0.2905551731698784, 0.2542357765236439, 0.0018500738518652474, 0.0018500738518652468, 0.021378633826174363, 0.17102907060939507, 0.14965043678322087, 0.5914157987683232, 2.492221334007053, 34.89109867609886, 69.78219735219757, 0.06838278542354957, 0.21882504381453613, 3.0635506134035166, 6.127101226807019, 0.020113671962140332, 0.05779012387438115, 0.8090617342413388, 1.6181234684826735, 0.008689838125755772, 0.02360203012214284, 0.3304284217100008, 0.66085684342, 0.004610860536255201, 0.012106465548744895, 0.16949051768242912, 0.3389810353648574, 0.002775110777797868, 0.007126211275391439, 0.09976695785548051, 0.19953391571096052, 224.9999999372875, 449.99999987457494, 224.9999999372875, 24.86291506178918, 49.72583012357836, 24.86291506178918, 3.9526280922242063, 7.905256184448414, 3.9526280922242063, 1.3819218057633562, 2.763843611526713, 1.3819218057633562, 0.6644090447915999, 1.3288180895831996, 0.6644090447915999, 0.378880959839935, 0.7577619196798702, 0.378880959839935, 236.24999995916613, 236.24999995916636, 39.37499999319434, 9.147859375871146, 45.99831263163245, 45.99831263163249, 7.666385438605405, 1.3393841253650922, 7.211194944079919, 7.2111949440799235, 1.2018658240133193, 0.4472278787683122, 2.494299377173945, 2.4942993771739466, 0.41571656286232406, 0.2088868497128045, 1.1898422770550403, 1.1898422770550408, 0.19830704617583989, 0.11680620884347706, 0.6745778866710682, 0.6745778866710688, 0.11242964777851132, 118.12499997958314, 314.9999999455548, 118.12499997958314, 18.295718751742296, 22.99915631581624, 61.331083508843264, 22.99915631581624, 2.6787682507301853, 3.6055974720399617, 9.614926592106558, 3.6055974720399617, 0.8944557575366244, 1.2471496885869733, 3.325732502898593, 1.2471496885869733, 0.41777369942560905, 0.5949211385275205, 1.5864563694067197, 0.5949211385275205, 0.23361241768695415, 0.33728894333553433, 0.8994371822280909, 0.33728894333553433, 39.37499999319434, 236.24999995916636, 236.24999995916613, 9.147859375871146, 7.666385438605405, 45.99831263163249, 45.99831263163245, 1.3393841253650922, 1.2018658240133193, 7.2111949440799235, 7.211194944079919, 0.4472278787683122, 0.41571656286232406, 2.4942993771739466, 2.494299377173945, 0.2088868497128045, 0.19830704617583989, 1.1898422770550408, 1.1898422770550403, 0.11680620884347706, 0.11242964777851132, 0.6745778866710688, 0.6745778866710682, 192.85714284335833, 128.5714285622389, 12.857142856223872, 5.8904120173888, 74.0147550944144, 49.343170062942946, 4.934317006294287, 0.791178878178647, 10.761063673481932, 7.174042448987955, 0.7174042448987946, 0.2512305823946569, 3.5525422613613777, 2.368361507574253, 0.23683615075742495, 0.1136000897091472, 1.643003148136361, 1.0953354320909074, 0.10953354320909058, 0.06212451421999867, 0.9115567078672215, 0.6077044719114811, 0.06077044719114802, 128.5714285622389, 205.71428569958246, 38.57142856867167, 5.890412017388805, 2.9452060086944014, 49.343170062942946, 78.94907210070878, 14.802951018882881, 0.7911788781786475, 0.3955894390893237, 7.174042448987955, 11.47846791838074, 2.1522127346963864, 0.251230582394657, 0.1256152911973285, 2.368361507574253, 3.7893784121188085, 0.7105084522722757, 0.11360008970914727, 0.056800044854573635, 1.0953354320909074, 1.7525366913454539, 0.32860062962727216, 0.062124514219998726, 0.031062257109999353, 0.6077044719114811, 0.9723271550583709, 0.18231134157344428, 77.14285713734337, 231.4285714120298, 77.1428571373434, 0.9817353362314665, 7.853882689851732, 0.9817353362314665, 29.605902037765773, 88.8177061132972, 29.605902037765784, 0.13186314636310786, 1.0549051709048627, 0.13186314636310786, 4.304425469392775, 12.913276408178309, 4.304425469392776, 0.0418717637324428, 0.33497410985954246, 0.0418717637324428, 1.4210169045445522, 4.2630507136336515, 1.4210169045445527, 0.018933348284857866, 0.1514667862788629, 0.018933348284857866, 0.6572012592545448, 1.9716037777636315, 0.6572012592545449, 0.010354085703333113, 0.08283268562666489, 0.010354085703333113, 0.3646226831468888, 1.0938680494406652, 0.3646226831468889, 38.57142856867168, 205.71428569958246, 128.57142856223894, 2.9452060086944014, 5.890412017388805, 14.802951018882887, 78.94907210070878, 49.34317006294296, 0.3955894390893237, 0.7911788781786475, 2.1522127346963873, 11.47846791838074, 7.174042448987958, 0.1256152911973285, 0.251230582394657, 0.710508452272276, 3.7893784121188085, 2.3683615075742535, 0.056800044854573635, 0.11360008970914727, 0.3286006296272723, 1.7525366913454539, 1.095335432090908, 0.031062257109999353, 0.062124514219998726, 0.1823113415734444, 0.9723271550583709, 0.6077044719114812, 12.857142856223872, 128.5714285622389, 192.85714284335822, 5.8904120173888, 4.934317006294287, 49.343170062942946, 74.01475509441437, 0.791178878178647, 0.7174042448987946, 7.174042448987955, 10.761063673481926, 0.2512305823946569, 0.23683615075742495, 2.368361507574253, 3.5525422613613773, 0.1136000897091472, 0.10953354320909058, 1.0953354320909074, 1.6430031481363603, 0.06212451421999867, 0.06077044719114802, 0.6077044719114811, 0.9115567078672211, 112.4999732620907, 56.249986631045466, 4.017856187931811, 2.856591823639904, 113.15155375264474, 56.57577687632249, 4.041126919737314, 0.3393237208189685, 13.619360601714082, 6.809680300857055, 0.48640573577550295, 0.10029975434198493, 4.044301419835181, 2.022150709917595, 0.14443933642268503, 0.043309780404644296, 1.7498668090899654, 0.8749334045449845, 0.06249524318178448, 0.02294641822319993, 0.9280551370437763, 0.4640275685218892, 0.03314482632299202, 84.37497994656822, 96.42854851036358, 12.053568563795437, 1.9043945490932699, 1.9043945490932699, 84.86366531448378, 96.98704607369562, 12.123380759211944, 0.22621581387931236, 0.22621581387931236, 10.214520451285587, 11.673737658612087, 1.4592172073265093, 0.06686650289465665, 0.06686650289465665, 3.033226064876393, 3.4665440741444455, 0.43331800926805525, 0.02887318693642954, 0.02887318693642954, 1.3124001068174775, 1.4998858363628296, 0.1874857295453535, 0.01529761214879996, 0.01529761214879996, 0.696041352782834, 0.7954758317518094, 0.09943447896897609, 60.2678428189772, 120.53568563795451, 24.10713712759087, 0.1904394549093268, 3.0470312785492353, 1.1426367294559625, 60.61690379605974, 121.23380759211959, 24.24676151842388, 0.02262158138793121, 0.3619453022069002, 0.13572948832758747, 7.29608603663255, 14.592172073265113, 2.9184344146530186, 0.006686650289465657, 0.10698640463145073, 0.04011990173679399, 2.166590046340277, 4.333180092680559, 0.8666360185361103, 0.002887318693642951, 0.04619709909828732, 0.017323912161857736, 0.9374286477267677, 1.8748572954535372, 0.3749714590907069, 0.0015297612148799944, 0.024476179438079962, 0.00917856728927998, 0.4971723948448805, 0.994344789689762, 0.19886895793795215, 40.178561879318124, 128.57139801381774, 40.17856187931816, 0.5713183647279808, 3.4279101883678824, 0.5713183647279811, 40.41126919737315, 129.31606143159382, 40.41126919737318, 0.0678647441637937, 0.4071884649827619, 0.06786474416379373, 4.864057357755032, 15.564983544816071, 4.864057357755035, 0.020059950868396993, 0.12035970521038182, 0.020059950868397, 1.444393364226851, 4.622058765525914, 1.444393364226852, 0.008661956080928861, 0.05197173648557313, 0.008661956080928866, 0.6249524318178451, 1.9998477818171, 0.6249524318178454, 0.004589283644639988, 0.0275357018678399, 0.004589283644639989, 0.3314482632299203, 1.0606344423357428, 0.33144826322992055, 24.107137127590878, 120.5356856379546, 60.2678428189772, 1.1426367294559625, 3.0470312785492353, 0.19043945490932673, 24.24676151842389, 121.23380759211966, 60.61690379605974, 0.13572948832758747, 0.3619453022069002, 0.022621581387931203, 2.9184344146530194, 14.592172073265122, 7.29608603663255, 0.04011990173679399, 0.10698640463145073, 0.006686650289465655, 0.8666360185361107, 4.333180092680561, 2.166590046340277, 0.017323912161857736, 0.04619709909828732, 0.00288731869364295, 0.374971459090707, 1.8748572954535385, 0.9374286477267677, 0.00917856728927998, 0.024476179438079962, 0.001529761214879994, 0.19886895793795217, 0.9943447896897627, 0.4971723948448805, 12.053568563795437, 96.42854851036358, 84.37497994656822, 1.9043945490932712, 1.9043945490932699, 12.123380759211944, 96.98704607369562, 84.86366531448378, 0.22621581387931247, 0.22621581387931236, 1.4592172073265093, 11.673737658612087, 10.214520451285587, 0.06686650289465668, 0.06686650289465665, 0.43331800926805525, 3.4665440741444455, 3.033226064876393, 0.02887318693642956, 0.02887318693642954, 0.1874857295453535, 1.4998858363628296, 1.3124001068174775, 0.015297612148799968, 0.01529761214879996, 0.09943447896897609, 0.7954758317518094, 0.696041352782834, 4.0178561879318035, 56.24998663104544, 112.49997326209059, 2.856591823639903, 4.041126919737305, 56.57577687632247, 113.15155375264466, 0.3393237208189684, 0.48640573577550206, 6.809680300857053, 13.619360601714073, 0.1002997543419849, 0.14443933642268478, 2.0221507099175935, 4.044301419835177, 0.04330978040464428, 0.06249524318178436,
0.874933404544984, 1.7498668090899636, 0.02294641822319992, 0.033144826322991955, 0.4640275685218889, 0.9280551370437756, 0.9052124319896028, 168.32462537140273, 67.32985014856123, 3.74054723047562, 0.08987860988976477, 12.63917944724771, 5.055671778899094, 0.2808706543832827, 0.02396619311529677, 3.016728922173512, 1.206691568869407, 0.06703842049274478, 0.009691570019648795, 1.150480613118033, 0.460192245247214, 0.025566235847067426, 0.004910344640351765, 0.5624576588157804, 0.22498306352631262, 0.012499059084795133, 0.4526062159948021, 0.6789093239922035, 134.65970029712267, 119.6975113752203, 11.221641691426893, 0.04493930494488248, 0.06740895741732374, 10.111343557798202, 8.98786094026508, 0.8426119631498507, 0.011983096557648407, 0.017974644836472615, 2.4133831377388177, 2.145229455767841, 0.20111526147823494, 0.004845785009824407, 0.007268677514736613, 0.9203844904944295, 0.8181195471061606, 0.0766987075412025, 0.002455172320175888, 0.0036827584802638326, 0.4499661270526258, 0.39996989071344574, 0.037497177254385505, 0.032329015428200096, 0.7758963702768031, 0.4849352314230016, 104.7353224533176, 157.1029836799768, 22.443283382853757, 0.0032099503532058847, 0.07703880847694132, 0.04814925529808829, 7.864378322731933, 11.796567484097928, 1.6852239262996984, 0.000855935468403456, 0.020542451241682972, 0.012839032026051847, 1.8770757737968575, 2.8156136606952935, 0.40223052295646927, 0.00034612750070174266, 0.008307060016841835, 0.005191912510526142, 0.7158546037178894, 1.0737819055768365, 0.15339741508240473, 0.00017536945144113446, 0.004208866834587232, 0.0026305417716170187, 0.3499736543742644, 0.524960481561398,
0.0749943545087709, 0.09698704628460027, 0.9698704628460042, 0.3232901542820012, 78.5514918399883, 179.5462670628304, 37.40547230475628, 0.009629851059617655, 0.09629851059617667, 0.03209950353205887, 5.898283742048957, 13.481791410397614, 2.8087065438328334, 0.0025678064052103676, 0.025678064052103718, 0.008559354684034567, 1.407806830347645, 3.2178441836517595, 0.6703842049274492, 0.0010383825021052276, 0.010383825021052295, 0.003461275007017429, 0.5368909527884176, 1.2271793206592405, 0.25566235847067476, 0.0005261083543234034, 0.005261083543234042, 0.0017536945144113459, 0.26248024078069865, 0.5999548360701682, 0.12499059084795161, 0.1939740925692006, 1.0345284937024015, 0.19397409256920062, 56.108208457134424, 187.0273615237809, 56.10820845713441, 0.019259702119235306, 0.10271841130258813, 0.019259702119235313, 4.21305981574925, 14.04353271916413, 4.21305981574925, 0.005135612810420737, 0.02738993498891055, 0.005135612810420739, 1.0055763073911739, 3.351921024637237, 1.0055763073911734, 0.002076765004210456, 0.011076080022455746, 0.0020767650042104566, 0.3834935377060122, 1.2783117923533704, 0.383493537706012,
0.0010522167086468067, 0.005611822446116295, 0.0010522167086468072, 0.18748588627192742, 0.6249529542397563, 0.18748588627192736, 0.3232901542820012, 0.9698704628460046, 0.0969870462846003, 37.40547230475628, 179.5462670628304, 78.55149183998836, 0.03209950353205888, 0.0962985105961767, 0.009629851059617656, 2.8087065438328334, 13.481791410397614, 5.898283742048961, 0.00855935468403457, 0.025678064052103728,
0.0025678064052103685, 0.6703842049274492, 3.2178441836517595, 1.407806830347646, 0.00346127500701743,
0.010383825021052298, 0.001038382502105228, 0.25566235847067476, 1.2271793206592405, 0.5368909527884181, 0.0017536945144113465, 0.005261083543234043, 0.0005261083543234036, 0.12499059084795161, 0.5999548360701682, 0.26248024078069887, 0.4849352314230016, 0.7758963702768031, 0.03232901542820006, 22.443283382853732, 157.1029836799768, 104.73532245331756, 0.04814925529808829, 0.07703880847694132, 0.0032099503532058817, 1.6852239262996966, 11.796567484097928, 7.86437832273193, 0.012839032026051847, 0.020542451241682972, 0.0008559354684034552, 0.40223052295646883, 2.8156136606952935, 1.8770757737968575, 0.005191912510526142, 0.008307060016841835, 0.0003461275007017422, 0.1533974150824046, 1.0737819055768365, 0.7158546037178892, 0.0026305417716170187, 0.004208866834587232, 0.00017536945144113432, 0.07499435450877082, 0.524960481561398, 0.3499736543742643, 0.6789093239922035, 0.4526062159948021, 11.221641691426903, 119.6975113752203, 134.6597002971226, 0.06740895741732374, 0.04493930494488246, 0.8426119631498513, 8.98786094026508, 10.111343557798198, 0.017974644836472615, 0.011983096557648404, 0.20111526147823508, 2.145229455767841, 2.413383137738817, 0.007268677514736613, 0.004845785009824405, 0.07669870754120256, 0.8181195471061606, 0.920384490494429, 0.0036827584802638326, 0.002455172320175887, 0.03749717725438553, 0.39996989071344574, 0.4499661270526256, 0.9052124319896017, 3.7405472304756184, 67.32985014856115, 168.32462537140268, 0.08987860988976468, 0.28087065438328257, 5.055671778899089, 12.639179447247702, 0.02396619311529674, 0.06703842049274472, 1.2066915688694058, 3.016728922173511, 0.009691570019648784, 0.02556623584706741, 0.4601922452472136, 1.1504806131180325, 0.00491034464035176, 0.012499059084795122, 0.22498306352631234, 0.56245765881578, 472.4999611813844, 944.9999223627688, 472.4999611813844, 45.698615732273254, 91.39723146454652, 45.698615732273254, 7.057053416997079, 14.114106833994159, 7.057053416997079, 2.4240860649990537, 4.848172129998107, 2.4240860649990537, 1.1527235812188221, 2.3054471624376447, 1.1527235812188221, 518.3999743128963, 518.3999743128967, 86.39999571881603, 19.685650812195185, 79.55089083649675, 79.5508908364968, 13.258481806082784, 2.8121865089841003, 12.278043429119872, 12.278043429119881, 2.046340571519978, 0.9231699170739912, 4.205860372901326, 4.205860372901328, 0.7009767288168872, 0.4262551558055232, 1.9944174236461918, 1.9944174236461931, 0.3324029039410319, 259.19998715644834, 691.1999657505282, 259.19998715644834, 39.371301624390384, 39.775445418248395, 106.06785444866233, 39.775445418248395, 5.624373017968201, 6.139021714559941, 16.370724572159826, 6.139021714559941, 1.8463398341479826, 2.102930186450664, 5.6078138305350995, 2.102930186450664, 0.8525103116110462, 0.9972087118230964, 2.6592232315282556, 0.9972087118230964, 86.39999571881603, 518.3999743128967, 518.3999743128963, 19.685650812195185, 13.258481806082784, 79.5508908364968, 79.55089083649675, 2.8121865089841003, 2.046340571519978, 12.278043429119881, 12.278043429119872, 0.9231699170739912, 0.7009767288168872, 4.205860372901328,
4.205860372901326, 0.4262551558055232, 0.3324029039410319, 1.9944174236461931, 1.9944174236461918, 468.64284870565336, 312.42856580376906, 31.24285658037686, 14.296729495807206, 121.0967123825081, 80.73114158833876, 8.073114158833864, 1.898998165868031, 17.964643923331284, 11.97642928222086, 1.1976429282220842, 0.596494536379557, 5.9869140090740265, 3.9912760060493526, 0.39912760060493463, 0.26745144886915223, 2.7838351984412326, 1.8558901322941554, 0.18558901322941526, 312.42856580376906, 499.88570528603105,
93.7285697411307, 14.296729495807218, 7.148364747903608, 80.73114158833876, 129.16982654134213, 24.219342476501623, 1.8989981658680328, 0.9494990829340162, 11.97642928222086, 19.162286851553393, 3.5929287846662574, 0.5964945363795575, 0.2982472681897787, 3.9912760060493526, 6.38604160967897, 1.1973828018148056, 0.26745144886915245, 0.13372572443457623, 1.8558901322941554, 2.969424211670652, 0.5567670396882466, 187.45713948226145, 562.3714184467837, 187.45713948226157, 2.3827882493012016, 19.062305994409606, 2.3827882493012016, 48.438684953003275, 145.31605485900965, 48.43868495300329, 0.31649969431133856, 2.5319975544907085, 0.31649969431133856, 7.185857569332518, 21.557572707997526, 7.18585756933252, 0.09941575606325953, 0.7953260485060761, 0.09941575606325953, 2.394765603629612, 7.184296810888827, 2.394765603629613, 0.04457524147819204, 0.35660193182553623, 0.04457524147819204, 1.1135340793764938, 3.3406022381294767, 1.1135340793764943, 93.72856974113073, 499.88570528603105, 312.42856580376906, 7.148364747903608, 14.296729495807218, 24.219342476501634, 129.16982654134213, 80.73114158833879, 0.9494990829340162, 1.8989981658680328, 3.592928784666259, 19.162286851553393, 11.976429282220863, 0.2982472681897787, 0.5964945363795575, 1.197382801814806, 6.38604160967897, 3.9912760060493535, 0.13372572443457623, 0.26745144886915245, 0.5567670396882468, 2.969424211670652, 1.855890132294156, 31.24285658037686, 312.42856580376906, 468.6428487056532, 14.296729495807206, 8.073114158833864, 80.73114158833876, 121.09671238250806, 1.898998165868031, 1.1976429282220842, 11.97642928222086, 17.964643923331277, 0.596494536379557, 0.39912760060493463, 3.9912760060493526, 5.986914009074025, 0.26745144886915223, 0.18558901322941526, 1.8558901322941554, 2.7838351984412317, 359.9999995363062, 179.99999976815351, 12.857142840582366, 8.415435551328036, 175.69070419863502, 87.84535209931771, 6.2746680070941085, 1.0126901133849284, 23.5496127957637, 11.774806397881875, 0.8410575998487038, 0.2995206392851393, 7.35818734999443, 3.679093674997223, 0.26279240535694404, 0.1290155474879426, 3.2747575034249303, 1.637378751712469, 0.11695562512231894, 269.9999996522303, 308.57142817397715, 38.57142852174711, 5.610290367552025, 5.610290367552025, 131.76802814897658, 150.59203217025882, 18.824004021282335, 0.6751267422566192, 0.6751267422566192, 17.662209596822823, 20.18538239636891, 2.5231727995461117, 0.19968042619009296, 0.19968042619009296, 5.518640512495836, 6.307017728566663, 0.7883772160708321, 0.08601036499196178, 0.08601036499196178, 2.4560681275687037, 2.8069350029356586, 0.35086687536695693, 192.85714260873561, 385.7142852174716, 77.1428570434942, 0.5610290367552018, 8.97646458808325, 3.366174220531216, 94.12002010641169, 188.24004021282357, 37.648008042564655, 0.06751267422566183, 1.080202787610592, 0.4050760453539717, 12.615863997730564, 25.231727995461153, 5.046345599092222, 0.019968042619009273, 0.31948868190414903, 0.1198082557140558, 3.9418860803541618,
7.8837721607083315, 1.5767544321416642, 0.008601036499196168, 0.13761658398713897, 0.05160621899517709, 1.7543343768347852, 3.5086687536695744, 0.7017337507339139, 128.57142840582372, 411.4285708986351, 128.5714284058238, 1.6830871102656073, 10.098522661593636, 1.6830871102656078, 62.74668007094112, 200.78937622701116, 62.74668007094116, 0.2025380226769857, 1.2152281360619133, 0.20253802267698578, 8.41057599848704, 26.913843195158478, 8.410575998487047, 0.05990412785702787, 0.359424767142167, 0.059904127857027895, 2.6279240535694406, 8.409356971422193, 2.627924053569443, 0.025803109497588524, 0.15481865698553104, 0.02580310949758854, 1.16955625122319, 3.7425800039142008, 1.1695562512231907, 77.14285704349423, 385.7142852174719, 192.85714260873561, 3.366174220531216, 8.97646458808325, 0.5610290367552017, 37.64800804256467, 188.24004021282366, 94.12002010641169, 0.4050760453539717, 1.080202787610592, 0.06751267422566183, 5.0463455990922235, 25.23172799546117, 12.615863997730564, 0.1198082557140558, 0.31948868190414903, 0.019968042619009266, 1.5767544321416647, 7.883772160708336, 3.9418860803541618, 0.05160621899517709, 0.13761658398713897, 0.008601036499196165, 0.7017337507339141, 3.5086687536695766, 1.7543343768347852, 38.57142852174711, 308.57142817397715, 269.9999996522303, 5.610290367552027, 5.610290367552025, 18.824004021282335, 150.59203217025882, 131.76802814897658, 0.6751267422566196, 0.6751267422566192, 2.5231727995461117, 20.18538239636891, 17.662209596822823, 0.19968042619009305, 0.19968042619009296, 0.7883772160708321, 6.307017728566663, 5.518640512495836, 0.08601036499196181, 0.08601036499196178, 0.35086687536695693, 2.8069350029356586, 2.4560681275687037, 12.857142840582341, 179.9999997681534, 359.9999995363059, 8.41543555132803, 6.274668007094098, 87.84535209931765, 175.69070419863488, 1.012690113384928, 0.8410575998487021, 11.774806397881868, 23.549612795763682, 0.2995206392851393, 0.2627924053569434, 3.679093674997221, 7.3581873499944255, 0.12901554748794258, 0.11695562512231872, 1.6373787517124678, 3.274757503424927, 202.500000376786, 81.00000015071458, 4.500000008373026, 4.030977325153006, 248.47334490086823, 99.38933796034749, 5.521629886685965, 0.42400691653479294, 26.835719491688284, 10.734287796675336, 0.596349322037518, 0.11556321683393261, 7.399479085092924, 2.9597916340371757, 0.1644328685576207, 0.047146852080466484, 3.0377574225881694, 1.2151029690352706, 0.06750572050195938, 162.00000030142937, 144.0000002679374, 13.50000002511912, 2.015488662576507, 3.023232993864761, 198.77867592069526, 176.69215637395155, 16.56488966005795, 0.21200345826739686, 0.31800518740109535, 21.468575593350696, 19.083178305200647, 1.7890479661125593, 0.05778160841696641, 0.08667241262544965, 5.919583268074359, 5.261851793843882, 0.49329860567286354, 0.023573426040233287, 0.03536013906034994, 2.4302059380705443, 2.1601830560627087, 0.20251716150587876, 126.000000234445, 189.000000351668, 27.000000050238203, 0.14396347589832162, 3.4551234215597235, 2.1594521384748258, 154.60563682720738, 231.90845524081166, 33.129779320115844, 0.015143104161956888, 0.36343449988696575, 0.22714656242935347, 16.69778101705054, 25.046671525575874, 3.5780959322251134, 0.004127257744069022, 0.09905418585765664, 0.06190886616103536, 4.60412031961339, 6.9061804794201, 0.9865972113457255, 0.0016838161457309454, 0.04041158749754274, 0.025257242185964198, 1.8901601740548668, 2.835240261082307, 0.40503432301175685, 94.50000017583386, 216.000000401906, 45.00000008373037, 0.4318904276949649, 4.318904276949655, 1.4396347589832172, 115.95422762040567, 265.0382345609273, 55.21629886685979, 0.04542931248587067, 0.45429312485870743, 0.15143104161956902, 12.52333576278792, 28.62476745780096, 5.963493220375194, 0.012381773232207065, 0.12381773232207083, 0.041272577440690246, 3.453090239710046, 7.892777690765819, 1.6443286855762107, 0.0050514484371928375, 0.05051448437192844, 0.01683816145730947, 1.4176201305411518, 3.240274584094061, 0.6750572050195953, 67.50000012559555, 225.00000041865124, 67.50000012559555, 0.8637808553899299, 4.606831228746285, 0.8637808553899301, 82.8244483002897, 276.0814943342982, 82.82444830028966, 0.09085862497174134, 0.4845793331826197, 0.09085862497174137, 8.94523983056279, 29.817466101875887, 8.945239830562786, 0.024763546464414133, 0.1320722478102085,
0.024763546464414143, 2.4664930283643156, 8.22164342788103, 2.4664930283643147, 0.010102896874385673, 0.05388211666339017, 0.010102896874385677, 1.012585807529393, 3.3752860250979677, 1.012585807529393, 45.00000008373037, 216.000000401906, 94.50000017583396, 1.4396347589832177, 4.318904276949657, 0.43189042769496505, 55.21629886685979, 265.0382345609273, 115.95422762040579, 0.15143104161956905, 0.45429312485870754, 0.04542931248587068, 5.963493220375194, 28.62476745780096, 12.52333576278793, 0.04127257744069026, 0.12381773232207087, 0.012381773232207068, 1.6443286855762107, 7.892777690765819, 3.4530902397100482, 0.01683816145730947, 0.05051448437192846, 0.0050514484371928375, 0.6750572050195953, 3.240274584094061, 1.4176201305411529, 27.000000050238175, 189.000000351668, 126.000000234445, 2.1594521384748258, 3.4551234215597235, 0.1439634758983215, 33.12977932011581, 231.90845524081166, 154.60563682720735, 0.22714656242935347, 0.36343449988696575, 0.015143104161956874, 3.5780959322251094, 25.046671525575874, 16.697781017050534, 0.06190886616103536, 0.09905418585765664, 0.004127257744069019, 0.9865972113457244, 6.9061804794201, 4.604120319613388, 0.025257242185964198, 0.04041158749754274, 0.0016838161457309437, 0.4050343230117564, 2.835240261082307, 1.8901601740548661, 13.500000025119128, 144.0000002679374, 162.00000030142928, 3.023232993864761, 2.015488662576506, 16.564889660057958, 176.69215637395155, 198.77867592069518, 0.31800518740109535, 0.21200345826739683, 1.7890479661125605, 19.083178305200647, 21.46857559335069, 0.08667241262544965, 0.0577816084169664, 0.4932986056728639, 5.261851793843882, 5.919583268074357, 0.03536013906034994, 0.023573426040233277, 0.20251716150587892, 2.1601830560627087, 2.4302059380705425, 4.500000008373024, 81.00000015071448, 202.5000003767859, 4.030977325153002, 5.521629886685964, 99.38933796034738, 248.4733449008681, 0.4240069165347925, 0.5963493220375177, 10.734287796675325, 26.835719491688273, 0.11556321683393249, 0.16443286855762057, 2.9597916340371717, 7.39947908509292, 0.04714685208046643, 0.06750572050195935, 1.215102969035269, 3.0377574225881685, 1.281060074308129, 345.80434908935746, 115.26811636311928, 5.239459834687229, 0.11120733278828383, 22.64815787261164, 7.549385957537223, 0.34315390716078203, 0.027027318090292008, 4.917807638028433, 1.639269212676146, 0.07451223693982464, 0.010228139409804665, 1.7525523602172528, 0.5841841200724184, 0.026553823639655316, 0.5124240297232528, 1.0248480594465073, 288.17029090779863, 209.57839338748974, 15.718379504061721, 0.04448293311531365, 0.08896586623062742, 18.873464893843085, 13.726156286431324, 1.0294617214823485, 0.01081092723611683, 0.021621854472233693, 4.0981730316903695, 2.980489477592995, 0.22353671081947443, 0.004091255763921876, 0.008182511527843763, 1.4604603001810479, 1.0621529455862158, 0.07966147091896612, 0.02846800165129179, 0.9109760528413408, 0.7971040462361716, 235.77569256092622, 282.93083107311185, 31.43675900812345, 0.0024712740619618665, 0.07908076998278002, 0.06919567373493238, 15.441925822235255, 18.530310986682334, 2.0589234429646974, 0.0006006070686731565, 0.019219426197541085, 0.016816997922848416, 3.353050662292123, 4.023660794750553, 0.4470734216389489, 0.0002272919868845484, 0.007273343580305575, 0.006364175632767367, 1.194922063784494, 1.433906476541395, 0.1593229418379323, 0.08540400495387566, 1.1956560693542608, 0.5978280346771297, 188.62055404874133, 335.32542941998526, 52.39459834687272, 0.007413822185885625, 0.10379351060239887, 0.05189675530119939, 12.35354065778823, 21.961850058290228, 3.4315390716078493, 0.0018018212060194757, 0.025225496884272697, 0.012612748442136333, 2.682440529833703, 4.7687831641488145, 0.7451223693982526, 0.0006818759606536474, 0.009546263449151078, 0.004773131724575533, 0.9559376510275972, 1.6994447129379535, 0.26553823639655544, 0.17080800990775102, 1.366464079262011, 0.42702002476937795, 146.7048753712432, 366.7621884281092, 78.59189752030886, 0.014827644371771222, 0.11862115497417003, 0.03706911092942808, 9.608309400501952, 24.02077350125496, 5.1473086074117615, 0.0036036424120389453, 0.02882913929631162, 0.009009106030097368, 2.0863426343151015, 5.215856585787773, 1.1176835540973762, 0.0013637519213072927, 0.010910015370458364, 0.0034093798032682333, 0.743507061910353, 1.8587676547758887, 0.3983073545948321, 0.28468001651291863, 1.423400082564589, 0.28468001651291863, 110.02865652843265, 377.2411080974812, 110.02865652843265, 0.024712740619618718, 0.12356370309809328, 0.024712740619618718, 7.206232050376481, 24.707081315576364, 7.206232050376481, 0.006006070686731579, 0.030030353433657814, 0.006006070686731579, 1.5647569757363298, 5.364881059667386, 1.5647569757363298, 0.002272919868845489, 0.011364599344227413, 0.002272919868845489, 0.5576302964327661, 1.9118753020551866, 0.5576302964327661, 0.4270200247693778, 1.3664640792620106, 0.17080800990775083, 78.59189752030883, 366.7621884281092, 146.70487537124308, 0.03706911092942807, 0.11862115497417, 0.014827644371771205, 5.14730860741176, 24.02077350125496, 9.608309400501945, 0.009009106030097367, 0.02882913929631161, 0.003603642412038941, 1.1176835540973757, 5.215856585787773, 2.0863426343151, 0.0034093798032682325, 0.010910015370458359, 0.0013637519213072912, 0.39830735459483196, 1.8587676547758887, 0.7435070619103523, 0.5978280346771301, 1.1956560693542608, 0.08540400495387569, 52.39459834687273, 335.32542941998537, 188.62055404874144, 0.051896755301199415, 0.10379351060239887, 0.007413822185885627, 3.43153907160785, 21.961850058290235, 12.353540657788235, 0.01261274844213634, 0.025225496884272697, 0.0018018212060194763, 0.7451223693982527, 4.768783164148816, 2.682440529833704, 0.0047731317245755375, 0.009546263449151078, 0.0006818759606536477, 0.26553823639655544, 1.699444712937954, 0.9559376510275976, 0.7971040462361716, 0.9109760528413408, 0.028468001651291774, 31.436759008123456, 282.93083107311173, 235.7756925609261, 0.06919567373493238, 0.07908076998278002, 0.0024712740619618648, 2.0589234429646983, 18.530310986682327, 15.441925822235248, 0.016816997922848416, 0.019219426197541085, 0.0006006070686731562, 0.4470734216389491, 4.023660794750551, 3.3530506622921212, 0.006364175632767367, 0.007273343580305575, 0.00022729198688454825, 0.15932294183793236, 1.4339064765413945, 1.1949220637844935, 1.0248480594465073, 0.5124240297232523, 15.718379504061728, 209.57839338748965, 288.17029090779863, 0.08896586623062742, 0.04448293311531359, 1.029461721482349, 13.72615628643132, 18.873464893843085, 0.021621854472233693, 0.01081092723611682, 0.22353671081947454, 2.9804894775929935, 4.0981730316903695, 0.008182511527843763, 0.004091255763921872, 0.07966147091896619, 1.0621529455862153, 1.4604603001810479, 1.281060074308129, 5.239459834687229, 115.26811636311928, 345.80434908935746, 0.11120733278828383, 0.34315390716078203, 7.549385957537223, 22.64815787261164, 0.027027318090292008, 0.07451223693982464, 1.639269212676146, 4.917807638028433, 0.010228139409804665, 0.026553823639655316, 0.5841841200724184, 1.7525523602172528, 881.9999999449774, 1763.9999998899548, 881.9999999449774, 77.42240967209004, 154.84481934418008, 77.42240967209004, 11.667440253688653, 23.33488050737731, 11.667440253688653, 3.944612895199086, 7.8892257903981715, 3.944612895199086, 992.2499999512207, 992.2499999512215, 165.37499999187006, 37.391581816442795, 128.76081568551652, 128.76081568551663, 21.46013594758608, 5.237749002797667, 19.53512359478478, 19.535123594784796, 3.2558539324641287, 1.6947808455658628, 6.617933761749877, 6.617933761749884, 1.1029889602916458, 496.12499997561065, 1322.9999999349611, 496.12499997561065, 74.78316363288556, 64.38040784275832, 171.6810875806887, 64.38040784275832, 10.475498005595332, 9.767561797392398, 26.046831459713047, 9.767561797392398, 3.3895616911317252, 3.308966880874941, 8.823911682333172, 3.308966880874941, 165.37499999187006, 992.2499999512215, 992.2499999512207, 37.391581816442795, 21.46013594758608, 128.76081568551663, 128.76081568551652, 5.237749002797667, 3.2558539324641287, 19.535123594784796, 19.53512359478478, 1.6947808455658628, 1.1029889602916458, 6.617933761749884, 6.617933761749877, 944.9999999698106, 629.9999999798739, 62.9999999979873, 29.43180821678096, 187.99696041618657, 125.33130694412444, 12.533130694412424, 3.871479873871827, 27.944274237104164, 18.629516158069446, 1.8629516158069415, 1.204477060405814, 9.31532095012401, 6.2102139667493415, 0.6210213966749333, 629.9999999798739, 1007.9999999677993, 188.99999999396215, 29.431808216780983, 14.715904108390488, 125.33130694412444, 200.53009111059927, 37.59939208323731, 3.8714798738718303, 1.9357399369359147, 18.629516158069446, 29.80722585291114, 5.588854847420832, 1.204477060405815, 0.6022385302029073, 6.2102139667493415, 9.936342346798957, 1.8630641900248024, 377.99999998792447, 1133.9999999637719, 377.9999999879246, 4.905301369463493, 39.242410955707946, 4.905301369463493, 75.19878416647467, 225.59635249942374, 75.1987841664747, 0.6452466456453045, 5.161973165162435, 0.6452466456453045, 11.177709694841672, 33.53312908452497, 11.177709694841676, 0.20074617673430234, 1.6059694138744187, 0.20074617673430234, 3.7261283800496066, 11.178385140148805, 3.7261283800496074, 188.99999999396223, 1007.9999999677993, 629.9999999798741, 14.715904108390488, 29.431808216780983, 37.59939208323733, 200.53009111059927, 125.33130694412444, 1.9357399369359147, 3.8714798738718303, 5.588854847420835, 29.80722585291114, 18.629516158069446, 0.6022385302029073, 1.204477060405815, 1.8630641900248028, 9.936342346798957, 6.210213966749343, 62.9999999979873, 629.9999999798739, 944.9999999698103, 29.43180821678096, 12.533130694412424, 125.33130694412444, 187.9969604161865, 3.871479873871827, 1.8629516158069415, 18.629516158069446, 27.944274237104153, 1.204477060405814, 0.6210213966749333, 6.2102139667493415, 9.315320950124006, 808.4999999868883, 404.2499999934451, 28.874999999531735, 19.429089813850627, 262.25094569890877, 131.12547284945467, 9.366105203532458, 2.3535912126274674, 36.7283767136629, 18.36418835683149, 1.3117277397736749, 0.6952531010198848, 11.746830503609779, 5.8734152518049, 0.41952966084320636, 606.3749999901678, 692.9999999887624, 86.62499999859523, 12.95272654256709, 12.95272654256709, 196.6882092741821, 224.78652488477928, 28.09831561059738, 1.5690608084183117, 1.5690608084183117, 27.54628253524724, 31.48146575456824, 3.935183219321026,
0.46350206734659, 0.46350206734659, 8.810122877707355, 10.068711860236967, 1.2585889825296197, 433.12499999297614, 866.2499999859532, 173.24999999719043, 1.2952726542567075, 20.72436246810736, 7.771635925540255, 140.49157805298694, 280.98315610597416, 56.19663122119476, 0.156906080841831, 2.5104972934693017, 0.9414364850509873, 19.67591609660514, 39.351832193210306, 7.87036643864205, 0.04635020673465895, 0.7416033077545447, 0.2781012404079541, 6.2929449126481, 12.585889825296212, 2.5171779650592385, 288.7499999953175, 923.999999985014, 288.74999999531764, 3.8858179627701253, 23.31490777662074, 3.885817962770127, 93.66105203532463, 299.7153665130382, 93.6610520353247, 0.47071824252549344, 2.824309455152959, 0.4707182425254936, 13.117277397736755, 41.97528767275754, 13.117277397736764, 0.139050620203977, 0.8343037212238613, 0.13905062020397704, 4.195296608432066, 13.424949146982582, 4.195296608432068, 173.24999999719046, 866.2499999859539, 433.12499999297614, 7.771635925540255, 20.72436246810736, 1.295272654256707, 56.19663122119477, 280.9831561059744, 140.49157805298694, 0.9414364850509873, 2.5104972934693017, 0.15690608084183097, 7.870366438642055, 39.35183219321034, 19.67591609660514, 0.2781012404079541, 0.7416033077545447, 0.04635020673465894, 2.5171779650592394, 12.585889825296219, 6.2929449126481, 86.62499999859523, 692.9999999887624, 606.3749999901678, 12.952726542567095, 12.95272654256709, 28.09831561059738, 224.78652488477928, 196.6882092741821, 1.5690608084183129, 1.5690608084183117, 3.935183219321026, 31.48146575456824, 27.54628253524724, 0.46350206734659033, 0.46350206734659, 1.2585889825296197, 10.068711860236967, 8.810122877707355, 28.874999999531678, 404.2499999934448, 808.4999999868877, 19.42908981385062, 9.366105203532438, 131.1254728494546, 262.25094569890854, 2.353591212627466, 1.3117277397736724, 18.364188356831477, 36.728376713662875, 0.6952531010198846, 0.41952966084320564, 5.873415251804898, 11.74683050360977, 601.3636363600964, 240.54545454403905, 13.363636363557708, 11.24637681879674, 357.1560586664259, 142.86242346657062, 7.936801303698359, 1.2236711369764697, 44.28292943569236, 17.713171774276983, 0.9840650985709423, 0.33781248964346927, 13.084860646070751, 5.233944258428313, 0.2907746810237947, 481.09090908807866, 427.6363636338483, 40.09090909067325, 5.623188409398381, 8.434782614097573, 285.7248469331417, 253.97764171834848, 23.81040391109515, 0.6118355684882361, 0.9177533527323545, 35.426343548554016, 31.49008315427027, 2.952195295712836, 0.1689062448217349, 0.25335936723260244, 10.467888516856638, 9.304789792761467, 0.872324043071387, 374.1818181796167, 561.2727272694264, 80.18181818134636, 0.4016563149570264, 9.639751558968646, 6.0248447243553995, 222.23043650355459, 333.3456547553327, 47.62080782219023, 0.04370254060630249, 1.048860974551261, 0.6555381090945377, 27.553822759986446, 41.330734139979775, 5.904390591425662, 0.012064731772981044, 0.2895535625515454, 0.18097097659471573, 8.141691068666272, 12.212536602999439, 1.744648086142771, 280.6363636347129, 641.4545454507721, 133.63636363557737, 1.2049689448710794, 12.04968944871081, 4.016563149570267, 166.6728273776661, 380.9664625775225, 79.36801303698377, 0.13110762181890745, 1.3110762181890765, 0.43702540606302515, 20.665367069989856, 47.235124731405385, 9.840650985709443, 0.03619419531894313, 0.36194195318943184, 0.12064731772981051, 6.106268301499711, 13.957184689142196, 2.907746810237954, 200.45454545336614, 668.1818181778851, 200.45454545336602, 2.409937889742159, 12.853002078624824, 2.4099378897421593, 119.05201955547567, 396.8400651849178, 119.05201955547562, 0.2622152436378149, 1.3984812994016773, 0.262215243637815, 14.760976478564169, 49.20325492854709, 14.760976478564166, 0.07238839063788627, 0.3860714167353927, 0.0723883906378863, 4.361620215356932, 14.538734051189733, 4.36162021535693, 133.63636363557737, 641.4545454507721, 280.636363634713, 4.016563149570269, 12.049689448710813, 1.2049689448710794, 79.36801303698377, 380.9664625775225, 166.67282737766624, 0.43702540606302526, 1.311076218189077, 0.1311076218189075, 9.840650985709443, 47.235124731405385, 20.665367069989873, 0.12064731772981055, 0.36194195318943195, 0.03619419531894314, 2.907746810237954, 13.957184689142196, 6.106268301499715, 80.1818181813463, 561.2727272694264, 374.1818181796166, 6.0248447243553995, 9.639751558968646, 0.401656314957026, 47.62080782219019, 333.3456547553327, 222.23043650355453, 0.6555381090945377, 1.048860974551261, 0.043702540606302444, 5.904390591425655, 41.330734139979775, 27.55382275998644, 0.18097097659471573, 0.2895535625515454, 0.01206473177298103, 1.7446480861427691, 12.212536602999439, 8.141691068666269, 40.09090909067328, 427.6363636338483, 481.09090908807855, 8.434782614097573, 5.623188409398379, 23.810403911095165, 253.97764171834848, 285.7248469331416, 0.9177533527323545, 0.6118355684882358, 2.9521952957128375, 31.49008315427027, 35.426343548553994, 0.25335936723260244, 0.16890624482173486, 0.8723240430713877, 9.304789792761467, 10.467888516856632, 13.363636363557703, 240.54545454403876, 601.3636363600962, 11.24637681879673, 7.9368013036983545, 142.86242346657045, 357.1560586664258, 1.2236711369764688, 0.9840650985709416, 17.71317177427696, 44.282929435692346, 0.3378124896434689, 0.2907746810237946, 5.2339442584283065, 13.084860646070748, 330.7501064559729, 110.25003548532446, 5.011365249332918, 5.3861719260738985, 479.1079559072378, 159.70265196907948, 7.25921145313996, 0.507752419795097, 46.86362460175436, 15.621208200584809, 0.7100549182083988, 0.1281100366980585, 12.030888758609063, 4.010296252869693, 0.18228619331225837, 275.62508871331147, 200.45460997331733, 15.03409574799879, 2.154468770429565, 4.308937540859136, 399.25662992269923, 290.3684581255993, 21.77763435941993, 0.20310096791803933, 0.4062019358360793, 39.05302050146208, 28.402196728336037, 2.130164754625201, 0.051244014679223536, 0.10248802935844718, 10.025740632174248, 7.291447732490355, 0.5468585799367763, 225.51143621998222, 270.6137234639791, 30.068191495997585, 0.11969270946830901, 3.830166702985903, 3.351395865112658, 326.6645153912996, 391.99741846956, 43.55526871883987, 0.011283387106557729, 0.3610683874098487, 0.31593483898361696, 31.95247131937807, 38.342965583253736, 4.260329509250402, 0.002846889704401304, 0.09110047054084205, 0.07971291172323665, 8.20287869905166, 9.843454438862004, 1.093717159873553, 180.4091489759861, 320.72737595730933, 50.113652493329596, 0.3590781284049283, 5.027093797669002, 2.5135468988344987, 261.3316123130401, 464.5895330009611, 72.59211453140018, 0.0338501613196733, 0.4739022584754269, 0.23695112923771314, 25.56197705550251, 45.44351476533787, 7.100549182084047, 0.008540669113203942, 0.11956936758485531, 0.0597846837924276, 6.56230295924134, 11.666316371984626, 1.8228619331225986, 140.31822698132248, 350.7955674533074, 75.17047873999421, 0.7181562568098553, 5.745250054478856, 1.7953906420246395, 203.25792068792, 508.14480171980176, 108.88817179710003, 0.06770032263934647, 0.541602581114773, 0.16925080659836633, 19.88153770983527, 49.70384427458836, 10.65082377312604, 0.01708133822640785, 0.1366507058112631, 0.04270334556601965, 5.104013412743262, 12.760033531858198, 2.7342928996838913, 105.2386702359921, 360.8182979519709, 105.2386702359921, 1.196927094683093, 5.984635473415449, 1.196927094683093, 152.44344051594035, 522.6632246260782, 152.44344051594035, 0.11283387106557755, 0.5641693553278861, 0.11283387106557755, 14.911153282376487, 51.12395411100482, 14.911153282376487, 0.02846889704401311, 0.14234448522006513, 0.02846889704401311, 3.828010059557455, 13.12460591848263, 3.828010059557455, 75.17047873999418, 350.7955674533074, 140.31822698132237, 1.7953906420246388, 5.745250054478853, 0.7181562568098545, 108.8881717971, 508.14480171980176, 203.2579206879198, 0.16925080659836628, 0.5416025811147728, 0.06770032263934642, 10.650823773126039, 49.70384427458836, 19.881537709835253, 0.04270334556601964, 0.13665070581126307, 0.017081338226407828, 2.73429289968389, 12.760033531858198, 5.1040134127432575, 50.11365249332961, 320.72737595730945, 180.4091489759862, 2.5135468988345, 5.027093797669002, 0.3590781284049285, 72.59211453140021, 464.58953300096124, 261.33161231304024, 0.2369511292377133, 0.4739022584754269, 0.03385016131967331, 7.100549182084047, 45.443514765337895, 25.561977055502517, 0.05978468379242764, 0.11956936758485531, 0.008540669113203943, 1.822861933122599, 11.66631637198463, 6.562302959241343, 30.068191495997596, 270.61372346397894, 225.5114362199821, 3.351395865112658, 3.830166702985903, 0.11969270946830894, 43.55526871883988, 391.9974184695598, 326.6645153912993, 0.31593483898361696, 0.3610683874098487, 0.011283387106557722, 4.260329509250404, 38.342965583253715, 31.95247131937806, 0.07971291172323665, 0.09110047054084205, 0.0028468897044013025, 1.0937171598735533, 9.843454438861999, 8.202878699051654, 15.034095747998798, 200.45460997331722, 275.62508871331147, 4.308937540859136,
2.1544687704295624, 21.777634359419938, 290.3684581255991, 399.25662992269923, 0.4062019358360793, 0.20310096791803908, 2.130164754625202, 28.40219672833602, 39.05302050146208, 0.10248802935844718, 0.051244014679223474, 0.5468585799367767, 7.2914477324903535, 10.025740632174248, 5.011365249332918, 110.25003548532446, 330.7501064559729, 5.3861719260738985, 7.25921145313996, 159.70265196907948, 479.1079559072378, 0.507752419795097, 0.7100549182083988, 15.621208200584809, 46.86362460175436, 0.1281100366980585, 0.18228619331225837, 4.010296252869693, 12.030888758609063, 1.7191583977373188, 635.9714897161454, 181.70613991889874, 6.98869768918838, 0.13242468337904043, 36.905104409274045, 10.544315545506876, 0.40555059790410874, 0.029510618707719073, 7.339104387439704, 2.096886967839917, 0.08064949876307337, 0.5730527992457739, 1.4326319981144369, 545.1184197566971, 335.45748908104207, 20.966093067565268, 0.04414156112634688, 0.11035390281586739, 31.63294663652068, 19.466428699397206, 1.2166517937123336, 0.009836872902573044, 0.024592182256432647, 6.290660903519762, 3.8711759406275186, 0.24194849628922155, 0.02604785451117148, 1.0419141804468621, 1.1721534530027213, 461.2540474864326, 461.25404748643354, 41.932186135130266, 0.0020064345966521262, 0.08025738386608529, 0.09028955684934603, 26.766339461671148, 26.766339461671205, 2.4333035874246516, 0.00044713058648059186, 0.017885223459223724, 0.020120876391626713, 5.322866918362835, 5.322866918362847, 0.48389699257843993, 0.07814356353351458, 1.4065841436032673, 0.9377227624021793, 384.37837290536174, 559.0958151350717, 69.8869768918838, 0.006019303789956389, 0.10834746821921537, 0.07223164547947701, 22.305282884726026, 32.44404783232877, 4.0555059790410874, 0.0013413917594417783, 0.02414505166995209, 0.016096701113301404, 4.435722431969045, 6.451959901045882, 0.8064949876307335,
0.1562871270670292, 1.6670626887149878, 0.7293399263128053, 314.4913960134778, 628.9827920269562, 104.83046533782607, 0.01203860757991278, 0.12841181418573708, 0.05618016870625984, 18.249776905684932, 36.4995538113699, 6.083258968561651, 0.0026827835188835566, 0.028616357534758097, 0.012519656421456637, 3.6292274443383077, 7.2584548886766225, 1.2097424814461044, 0.2604785451117169, 1.8233498157820198, 0.5470049447346054, 251.59311681078287, 670.9149781620894, 146.7626514729571, 0.020064345966521423, 0.14045042176565006, 0.04213512652969497, 14.599821524547982, 38.932857398794724, 8.516562555986345, 0.004471305864805955, 0.031299141053641705, 0.0093897423160925, 2.903381955470654, 7.742351881255098, 1.693639474024553, 0.3907178176675744, 1.8754455248043505, 0.39071781766757424, 195.68353529727565, 684.8923735404622, 195.68353529727565, 0.03009651894978206, 0.14446329095895338, 0.030096518949782054, 11.355416741315102, 39.74395859460271, 11.355416741315102, 0.006706958797208914, 0.03219340222660268, 0.006706958797208912, 2.2581859653660654, 7.9036508787812, 2.2581859653660654, 0.5470049447346054, 1.8233498157820198, 0.26047854511171714, 146.7626514729571, 670.9149781620895, 251.59311681078282, 0.04213512652969497, 0.14045042176565006, 0.02006434596652144, 8.516562555986345, 38.93285739879474, 14.599821524547979, 0.0093897423160925, 0.031299141053641705, 0.004471305864805957, 1.693639474024553, 7.742351881255101, 2.9033819554706533, 0.729339926312805, 1.6670626887149884, 0.1562871270670292, 104.83046533782611, 628.9827920269562, 314.4913960134778, 0.05618016870625981, 0.1284118141857371, 0.01203860757991278, 6.083258968561654, 36.4995538113699, 18.249776905684932, 0.012519656421456632, 0.028616357534758104, 0.0026827835188835566, 1.209742481446105, 7.2584548886766225, 3.6292274443383086, 0.9377227624021793, 1.4065841436032667, 0.07814356353351462, 69.88697689188378, 559.095815135072, 384.37837290536174, 0.07223164547947701, 0.10834746821921534, 0.006019303789956393, 4.055505979041086, 32.44404783232879, 22.305282884726026, 0.016096701113301404, 0.024145051669952074, 0.0013413917594417787, 0.8064949876307334, 6.451959901045885, 4.435722431969045, 1.172153453002721, 1.0419141804468617, 0.02604785451117148, 41.932186135130266, 461.25404748643336, 461.25404748643234, 0.09028955684934599, 0.08025738386608526, 0.0020064345966521262, 2.4333035874246516, 26.766339461671194, 26.76633946167113, 0.020120876391626706, 0.017885223459223717, 0.00044713058648059186, 0.48389699257843993, 5.322866918362845, 5.3228669183628305, 1.4326319981144369, 0.5730527992457739, 20.96609306756528, 335.45748908104196, 545.1184197566971, 0.11035390281586739, 0.04414156112634688, 1.2166517937123336, 19.466428699397206, 31.63294663652068, 0.024592182256432647, 0.009836872902573044, 0.2419484962892216, 3.8711759406275177, 6.290660903519762, 1.7191583977373193, 6.98869768918838, 181.70613991889874, 635.9714897161454, 0.13242468337904045, 0.40555059790410874, 10.544315545506876, 36.905104409274045, 0.029510618707719084, 0.08064949876307337, 2.096886967839917, 7.339104387439704, 1511.999981336902, 3023.9999626738045, 1511.999981336902, 123.28073323740263, 246.56146647480526,
123.28073323740263, 18.199190590484626, 36.39838118096926, 18.199190590484626, 1727.9999856380234, 1727.999985638025, 287.99999760633716, 64.9396366059751, 197.8669829415247, 197.86698294152484, 32.97783049025411, 8.952551403963264, 29.528692725943884, 29.528692725943902, 4.921448787657311, 863.9999928190123, 2303.9999808506977, 863.9999928190123, 129.8792732119502, 98.93349147076242, 263.8226439220329, 98.93349147076242, 17.90510280792653, 14.76434636297195, 39.37159030125851, 14.76434636297195, 287.99999760633716, 1727.999985638025, 1727.9999856380234, 64.9396366059751, 32.97783049025411, 197.86698294152484,
197.8669829415247, 8.952551403963264, 4.921448787657311, 29.528692725943902, 29.528692725943884, 1697.1428506283726, 1131.4285670855818, 113.14285670855801, 54.161101275743185, 279.62132817027776, 186.41421878018528, 18.641421878018498, 7.0647673423486586, 41.344065176591556, 27.562710117727715, 2.756271011772767, 1131.4285670855818, 1810.285707336933, 339.42857012567447, 54.16110127574323, 27.08055063787161, 186.41421878018528, 298.2627500482967, 55.924265634055565, 7.064767342348666, 3.532383671174332, 27.562710117727715, 44.10033618836438, 8.268813035318312, 678.8571402513493, 2036.5714207540454, 678.8571402513496, 9.026850212623865, 72.21480170099093, 9.026850212623865, 111.84853126811119, 335.54559380433307, 111.84853126811123, 1.1774612237247766, 9.419689789798213, 1.1774612237247766, 16.537626070636634, 49.61287821190984, 16.537626070636634, 339.4285701256746, 1810.285707336933, 1131.4285670855822, 27.08055063787161, 54.16110127574323, 55.92426563405558, 298.2627500482967, 186.4142187801853, 3.532383671174332, 7.064767342348666, 8.268813035318317, 44.10033618836438, 27.562710117727722, 113.14285670855801, 1131.4285670855818, 1697.142850628372, 54.161101275743185, 18.641421878018498, 186.41421878018528, 279.62132817027765, 7.0647673423486586, 2.756271011772767, 27.562710117727715, 41.344065176591535, 1535.9999989075973, 767.9999994538005, 54.8571428181285, 38.60125481886351, 378.2781385474956, 189.13906927374828, 13.509933519553421, 4.6920729554442815, 53.937108636483146, 26.968554318241633, 1.9263253084458276, 1151.9999991807013, 1316.5714276350855, 164.57142845438554, 25.734169879242344, 25.734169879242344, 283.7086039106225, 324.23840446928244, 40.52980055866026, 3.128048636962856, 3.128048636962856, 40.45283147736247, 46.23180740269992, 5.778975925337483, 822.857142271928, 1645.714284543857, 329.14285690877097, 2.5734169879242312, 41.17467180678779, 15.44050192754541, 202.64900279330138, 405.2980055866032, 81.05960111732051, 0.31280486369628524, 5.004877819140574, 1.8768291821777143, 28.89487962668742, 57.7897592533749, 11.557951850674964, 548.5714281812852, 1755.4285701801089, 548.5714281812856, 7.720250963772702,
46.32150578263617, 7.720250963772704, 135.09933519553425, 432.31787262570873, 135.09933519553434, 0.9384145910888565, 5.630487546533136, 0.9384145910888569, 19.263253084458277, 61.64240987026637, 19.263253084458295, 329.1428569087711, 1645.7142845438586, 822.857142271928, 15.44050192754541, 41.17467180678779, 2.5734169879242303, 81.05960111732055, 405.2980055866035, 202.64900279330138, 1.8768291821777143, 5.004877819140574, 0.31280486369628513, 11.557951850674971, 57.789759253374946, 28.89487962668742, 164.57142845438554, 1316.5714276350855, 1151.9999991807013, 25.73416987924235, 25.734169879242344, 40.52980055866026, 324.23840446928244, 283.7086039106225, 3.1280486369628577, 3.128048636962856, 5.778975925337483, 46.23180740269992, 40.45283147736247, 54.85714281812838, 767.9999994538, 1535.9999989075961, 38.60125481886348, 13.509933519553393, 189.13906927374816, 378.27813854749536, 4.6920729554442815, 1.9263253084458234, 26.968554318241623, 53.93710863648311, 1276.3636370206984, 510.54545480828045, 28.363636378237775, 24.979318050641833, 500.25892943760164, 200.10357177504108, 11.116865098613381, 2.7774596531218254, 66.13943941479066, 26.45577576591632, 1.469765320328683, 1021.0909096165625, 907.6363641036122, 85.09090913471357, 12.48965902532094, 18.734488537981413, 400.2071435500828, 355.73968315562956, 33.35059529584025, 1.3887298265609154, 2.0830947398413735, 52.911551531832714, 47.03249025051803, 4.409295960986062, 794.1818185906594, 1191.2727278859923, 170.1818182694269, 0.8921185018086367, 21.410844043407305, 13.381777527129557, 311.2722227611753, 466.9083341417641, 66.70119059168039, 0.09919498761149376, 2.380679702675853, 1.4879248141724069, 41.15342896920321, 61.73014345380497, 8.81859192197211, 595.6363639429953, 1361.4545461554176, 283.6363637823784, 2.67635550542591, 26.76355505425914, 8.921185018086375, 233.4541670708818, 533.6095247334441, 111.16865098613408, 0.29758496283448127, 2.975849628344817, 0.9919498761149383, 30.865071726902443, 70.54873537577701, 14.697653203286862, 425.4545456735677, 1418.181818911888, 425.4545456735675, 5.35271101085182, 28.547792057876332, 5.352711010851822, 166.75297647920112, 555.8432549306688, 166.75297647920104, 0.5951699256689625, 3.174239603567795, 0.5951699256689629, 22.046479804930293, 73.48826601643412, 22.046479804930286, 283.6363637823784, 1361.4545461554176, 595.6363639429957, 8.921185018086376, 26.763555054259147, 2.676355505425911, 111.16865098613408, 533.6095247334441,
233.45416707088197, 0.9919498761149386, 2.9758496283448177, 0.29758496283448127, 14.697653203286862, 70.54873537577701, 30.865071726902475, 170.1818182694267, 1191.2727278859923, 794.1818185906592, 13.381777527129557, 21.410844043407305, 0.8921185018086359, 66.7011905916803, 466.9083341417641, 311.2722227611753, 1.4879248141724069, 2.380679702675853, 0.09919498761149366, 8.8185919219721, 61.73014345380497, 41.1534289692032, 85.09090913471364, 907.6363641036122, 1021.090909616562, 18.734488537981413, 12.489659025320936, 33.350595295840264, 355.73968315562956, 400.2071435500826, 2.0830947398413735, 1.388729826560915, 4.409295960986065, 47.03249025051803, 52.911551531832686, 28.363636378237757, 510.54545480828, 1276.363637020698, 24.979318050641808, 11.116865098613376, 200.10357177504088, 500.2589294376015, 2.7774596531218223, 1.4697653203286818, 26.455775765916286, 66.13943941479063, 930.4615389236735, 310.15384630789157, 14.097902104904133, 14.421698650270484, 652.1820161603627, 217.39400538678782, 9.881545699399426, 1.4308512321103684, 75.09304011634393, 25.03101337211468, 1.1377733350961192, 775.38461576973, 563.9160841961668, 42.29370631471249, 5.768679460108208, 11.537358920216432, 543.4850134669704, 395.2618279759782, 29.644637098198345, 0.5723404928441486, 1.1446809856882991, 62.577533430286785, 45.5109334038449, 3.413320005288366, 634.4055947206886, 761.2867136648272, 84.58741262942499, 0.3204821922282334, 10.255430151303509, 8.973501382390552, 444.66955647297596, 533.6034677675717, 59.2892741963967, 0.03179669404689711, 1.0174942095007111, 0.8903074333131207, 51.19980007932557, 61.43976009519076, 6.826640010576732, 507.5244757765517, 902.2657347138716, 140.97902104904247, 0.9614465766847033, 13.460252073585867, 6.730126036792926, 355.7356451783814, 632.4189247615681, 98.81545699399507, 0.09539008214069165, 1.335461149969685, 0.6677305749848418, 40.95984006346054, 72.81749344615221, 11.377733350961286, 394.74125893731787, 986.853147343298, 211.46853157356318, 1.922893153369403, 15.38314522695526, 4.807232883423512, 276.6832795831854, 691.7081989579658, 148.2231854909922, 0.19078016428138295, 1.5262413142510671, 0.47695041070345784, 31.85765338269151, 79.64413345672905, 17.066600026441886, 296.0559442029891, 1015.0489515530996, 296.0559442029891, 3.204821922282341, 16.02410961141166, 3.204821922282341, 207.5124596873895, 711.4712903567602, 207.5124596873895, 0.31796694046897184, 1.589834702344855, 0.31796694046897184, 23.893240037018693, 81.91968012692075, 23.893240037018693, 211.4685315735631, 986.853147343298, 394.74125893731747, 4.80723288342351, 15.383145226955254, 1.9228931533694011, 148.22318549099217, 691.7081989579658, 276.68327958318514, 0.4769504107034577, 1.5262413142510665, 0.19078016428138275, 17.06660002644188, 79.64413345672905, 31.85765338269148, 140.97902104904253, 902.2657347138717, 507.52447577655204, 6.730126036792929, 13.460252073585867, 0.9614465766847038, 98.8154569939951, 632.4189247615683, 355.7356451783815, 0.6677305749848421, 1.335461149969685, 0.0953900821406917, 11.37773335096129, 72.81749344615224, 40.95984006346056, 84.58741262942502, 761.286713664827, 634.4055947206881, 8.973501382390552, 10.255430151303509, 0.32048219222823326, 59.28927419639671, 533.6034677675716, 444.6695564729756, 0.8903074333131207, 1.0174942095007111, 0.0317966940468971, 6.8266400105767335, 61.439760095190735, 51.19980007932554, 42.293706314712516, 563.9160841961667, 775.38461576973, 11.537358920216432, 5.768679460108202, 29.644637098198352, 395.261827975978, 543.4850134669704, 1.1446809856882991, 0.5723404928441479, 3.413320005288367, 45.51093340384488, 62.577533430286785, 14.097902104904133, 310.15384630789157, 930.4615389236735, 14.421698650270484, 9.881545699399426, 217.39400538678782, 652.1820161603627, 1.4308512321103684, 1.1377733350961192, 25.03101337211468, 75.09304011634393, 504.00000009225363, 144.00000002635826, 5.538461539475292, 6.925286706579552, 841.8052269417411, 240.51577912621195, 9.250606889469648, 0.5910953756661207, 75.18502158610283, 21.481434738886538, 0.8262090284187091, 432.0000000790754, 265.84615389481377, 16.615384618425978, 2.308428902193188, 5.77107225548298, 721.5473373786368, 444.02913069454263, 27.7518206684091, 0.19703179188870726, 0.4925794797217689, 64.4443042166597, 39.658033364098, 2.478627085256142, 365.5384616053688, 365.5384616053696, 33.230769236851735, 0.1049285864633265, 4.197143458533071, 4.72178639084971, 610.540054704996, 610.5400547049973, 55.503641336817836, 0.008955990540395763, 0.35823962161583156, 0.403019574317811, 54.52979587563473, 54.529795875634846, 4.95725417051225, 304.61538467114167, 443.07692315802433, 55.38461539475291, 0.3147857593899801, 5.66614366901966, 3.7774291126797777, 508.7833789208316, 740.0485511575733, 92.50606889469645, 0.026867971621187335, 0.4836234891813738, 0.32241565945424955, 45.4414965630291, 66.09672227349688, 8.262090284187092, 249.2307692763886, 498.46153855277765, 83.07692309212965, 0.6295715187799603, 6.715429533652947, 2.9380004209731574, 416.2773100261348, 832.5546200522706, 138.75910334204517, 0.0537359432423747, 0.5731833945853332, 0.2507677351310827, 37.17940627884198, 74.35881255768403, 12.393135426280677, 199.38461542111142, 531.6923077896319, 116.30769232898197, 1.0492858646332737, 7.345001052432918, 2.2035003157298734, 333.0218480209088, 888.0582613890924, 194.26274467886398, 0.08955990540395835, 0.626919337827709, 0.18807580134831248, 29.74352502307366, 79.31606672819665, 17.350389596793015, 155.07692310530894, 542.7692308685794, 155.07692310530894, 1.5739287969499063, 7.554858225359523, 1.5739287969499056, 259.01699290515137, 906.5594751680263, 259.01699290515137, 0.1343398581059372, 0.6448313189084962, 0.13433985810593715, 23.133852795723968, 80.96848478503358, 23.133852795723968, 116.30769232898197, 531.6923077896321, 199.38461542111136, 2.2035003157298734, 7.345001052432918, 1.0492858646332741, 194.26274467886398, 888.0582613890926, 333.0218480209087, 0.18807580134831248, 0.626919337827709, 0.08955990540395842, 17.350389596793015, 79.31606672819667, 29.743525023073648, 83.07692309212968, 498.46153855277765, 249.23076927638866, 2.9380004209731556, 6.715429533652949, 0.6295715187799603, 138.75910334204522, 832.5546200522706, 416.27731002613496, 0.2507677351310826, 0.5731833945853335, 0.0537359432423747, 12.393135426280683, 74.35881255768403, 37.17940627884199, 55.3846153947529, 443.0769231580245, 304.61538467114167, 3.7774291126797777, 5.666143669019657, 0.3147857593899802, 92.50606889469643, 740.0485511575735, 508.7833789208316, 0.32241565945424955, 0.48362348918137366, 0.026867971621187356, 8.262090284187089, 66.09672227349688, 45.4414965630291, 33.230769236851735, 365.5384616053695, 365.5384616053685, 4.721786390849708, 4.197143458533069, 0.1049285864633265, 55.503641336817836, 610.540054704997, 610.5400547049954, 0.4030195743178107, 0.35823962161583145, 0.008955990540395763, 4.95725417051225, 54.52979587563483, 54.52979587563469, 16.615384618425978, 265.84615389481377, 432.0000000790754, 5.77107225548298, 2.308428902193188, 27.75182066840911, 444.02913069454263, 721.5473373786368, 0.4925794797217689, 0.19703179188870726, 2.4786270852561425, 39.65803336409799, 64.4443042166597, 5.538461539475292, 144.00000002635826, 504.00000009225363, 6.925286706579555, 9.250606889469648, 240.51577912621195, 841.8052269417411, 0.5910953756661208, 0.8262090284187091, 21.481434738886538, 75.18502158610283, 2.2196035154623055, 1078.5759545108008, 269.6439886277012, 8.988132954256715, 0.15356418295207802, 56.160615844059066, 14.040153961014823, 0.46800513203382776, 0.6341724329892304, 1.902517298967694, 943.7539601969545, 503.335445438374, 26.964398862770267, 0.04387548084345087, 0.1316264425303528, 49.14053886355189, 26.20828739389425, 1.4040153961014894, 0.02439124742266257, 1.1707798762878034, 1.6098223298957284, 817.9200988373536, 701.0743704320182, 53.92879772553981, 0.00168751849397887, 0.08100088771098576, 0.11137622060260534, 42.588467015077946, 36.504400298638274, 2.808030792202942, 0.07317374226798823, 1.6098223298957324, 1.3415186082464454, 701.0743704320209, 862.8607636086382, 89.88132954256652, 0.005062555481936646, 0.11137622060260563, 0.09281351716883814, 36.504400298638416, 44.92849267524714, 4.6800513203382454, 0.14634748453597543, 1.9512997938130108, 1.0976061340198187, 593.2167749809383, 988.6946249682303, 134.8219943138496, 0.010125110963873224, 0.13500147951830999, 0.07593833222904936, 30.888338714232386, 51.48056452372062, 7.020076980507359, 0.24391247422662596, 2.1952122680396395, 0.8780849072158571, 494.34731248411543, 1078.5759545107985, 188.7507920393894, 0.01687518493978872, 0.1518766644580989, 0.060750665783239656, 25.74028226186033, 56.160615844058945, 9.828107772710299, 0.36586871133994014, 2.3415597525756255, 0.6829549278345559, 404.4659829415496, 1132.5047522363398, 251.6677227191868, 0.02531277740968316, 0.16200177542197283, 0.047250517831408635, 21.060230941522114, 58.96864663626197, 13.104143696947116, 0.5122161958759182, 2.390342247420937, 0.5122161958759182, 323.57278635323985, 1150.481018144855, 323.57278635323985, 0.035437888373556554, 0.1653768124099296, 0.035437888373556554, 16.84818475321771, 59.90465690032973, 16.84818475321771, 0.6829549278345559, 2.341559752575626, 0.36586871133994014, 251.66772271918694, 1132.5047522363407, 404.4659829415496, 0.047250517831408635, 0.1620017754219729, 0.025312777409683167, 13.104143696947126, 58.968646636262015, 21.060230941522114, 0.8780849072158567, 2.1952122680396395, 0.24391247422662585, 188.7507920393893, 1078.575954510798, 494.34731248411526, 0.060750665783239614, 0.1518766644580989, 0.016875184939788707, 9.828107772710295, 56.16061584405892, 25.74028226186031, 1.097606134019819, 1.9512997938130119, 0.14634748453597543, 134.8219943138496, 988.6946249682296, 593.2167749809383, 0.07593833222904942, 0.13500147951831007, 0.010125110963873224, 7.020076980507359, 51.48056452372059, 30.888338714232386, 1.3415186082464454, 1.6098223298957317, 0.07317374226798823, 89.88132954256642, 862.860763608638, 701.0743704320208, 0.09281351716883814, 0.1113762206026056, 0.005062555481936648, 4.680051320338242, 44.92849267524712, 36.504400298638416, 1.609822329895728, 1.170779876287803, 0.024391247422662576, 53.92879772553976, 701.074370432018, 817.9200988373536, 0.1113762206026053, 0.08100088771098576, 0.001687518493978871, 2.8080307922029397, 36.50440029863827, 42.588467015077946, 1.902517298967694, 0.6341724329892304, 26.964398862770274, 503.335445438374, 943.7539601969545, 0.1316264425303528, 0.04387548084345087, 1.4040153961014905, 26.20828739389425, 49.14053886355189, 2.2196035154623055, 8.988132954256715, 269.6439886277012, 1078.5759545108003, 0.15356418295207802, 0.46800513203382776, 14.040153961014823, 56.16061584405905, 2430.0000000166747, 4860.000000033349, 2430.0000000166747, 186.9429021059229, 373.8858042118458, 186.9429021059229, 2806.650000020477, 2806.650000020479, 467.7750000034125, 105.42677198079899, 291.6159537576239, 291.61595375762414, 48.60265895960397, 1403.3250000102391, 3742.200000027302, 1403.3250000102391, 210.853543961598, 145.80797687881207, 388.8212716768318, 145.80797687881207, 467.7750000034125, 2806.650000020479, 2806.650000020477, 105.42677198079899, 48.60265895960397, 291.61595375762414, 291.6159537576239, 2811.8571428765313, 1874.571428584355, 187.4571428584352, 91.85758685864805, 401.4133566667991, 267.6089044445328, 26.760890444453235, 1874.571428584355, 2999.3142857349712, 562.3714285753064, 91.85758685864812, 45.92879342932405, 267.6089044445328, 428.17424711125295, 80.28267133335983, 1124.7428571506134, 3374.2285714518357, 1124.7428571506139, 15.309597809774676, 122.4767824781974, 15.309597809774676, 160.56534266671972, 481.6960280001585, 160.56534266671977, 562.3714285753065, 2999.3142857349712, 1874.5714285843555, 45.92879342932405, 91.85758685864812, 80.28267133335986, 428.17424711125295, 267.6089044445328, 187.4571428584352, 1874.571428584355, 2811.8571428765304, 91.85758685864805, 26.760890444453235, 267.6089044445328, 401.41335666679896, 2632.500000013541, 1316.2500000067735, 94.01785714334078, 69.17874504897341, 529.7951002835817, 264.89755014179144, 18.921253581556492, 1974.375000010161, 2256.4285714401817, 282.05357143002243, 46.11916336598228, 46.11916336598228, 397.3463252126873, 454.1100859573565, 56.7637607446695, 1410.2678571501124, 2820.5357143002275, 564.1071428600447, 4.611916336598223, 73.79066138557172, 27.67149801958938, 283.81880372334757, 567.6376074466957, 113.52752148933897, 940.1785714334083, 3008.5714285869, 940.1785714334087, 13.835749009794682, 83.01449405876804, 13.83574900979469, 189.21253581556502, 605.480114609807, 189.21253581556516, 564.1071428600449, 2820.5357143002298, 1410.2678571501124, 27.67149801958938, 73.79066138557172, 4.611916336598221, 113.52752148933901, 567.637607446696, 283.81880372334757, 282.05357143002243, 2256.4285714401817, 1974.375000010161, 46.1191633659823, 46.11916336598228, 56.7637607446695, 454.1100859573565, 397.3463252126873, 94.01785714334058, 1316.250000006773, 2632.500000013539, 69.17874504897338, 18.921253581556456, 264.8975501417913, 529.7951002835814, 2319.5454545513708, 927.8181818205501, 51.54545454558606, 48.21888649633745, 684.2713750972988, 273.70855003892007, 15.206030557717765, 1855.636363641103, 1649.4545454587599, 154.63636363675863, 24.10944324816877, 36.16416487225317, 547.4171000778409, 486.5929778469703, 45.61809167315343, 1443.2727272764128, 2164.9090909146244, 309.2727272735168, 1.7221030891549083, 41.33047413971786, 25.831546337323637, 425.76885561609834, 638.653283424149, 91.23618334630673, 1082.4545454573106, 2474.181818188139, 515.4545454558618, 5.166309267464726, 51.66309267464733, 17.2210308915491, 319.3266417120742, 729.8894667704551, 152.06030557717799,
773.1818181837926, 2577.2727272793018, 773.1818181837924, 10.33261853492945, 55.10729885295698, 10.332618534929455, 228.09045836576698, 760.3015278858878, 228.09045836576692, 515.4545454558618, 2474.181818188139, 1082.4545454573117, 17.221030891549105, 51.66309267464735, 5.166309267464726, 152.06030557717799, 729.8894667704551, 319.32664171207443, 309.27272727351647, 2164.9090909146244, 1443.2727272764123, 25.831546337323637, 41.33047413971786, 1.7221030891549065, 91.23618334630663, 638.653283424149, 425.76885561609816, 154.63636363675874, 1649.4545454587599, 1855.6363636411018, 36.16416487225317, 24.10944324816876, 45.618091673153465, 486.5929778469703, 547.4171000778406, 51.54545454558603, 927.818181820549, 2319.5454545513694, 48.21888649633738, 15.206030557717755, 273.7085500389197, 684.2713750972983, 1892.5964267070576, 630.8654755690198, 28.675703434955388, 31.06322383922555, 871.8968941553063, 290.63229805176917, 13.210559002353111, 1577.163688922552, 1147.028137398219, 86.02711030486635, 12.425289535690247,
24.85057907138053, 726.5807451294239, 528.4223600941261, 39.63167700705943, 1290.4066545729977, 1548.4879854875992, 172.05422060973274, 0.6902938630939018, 22.08940361900494, 19.328228166629287, 594.4751551058924, 713.370186127072, 79.26335401411886, 1032.3253236584, 1835.2450198371594, 286.7570343495562, 2.070881589281713, 28.992342249944016, 14.496171124971992, 475.58012408471495, 845.4757761506058, 132.10559002353222, 802.9196961787552, 2007.2992404468948, 430.1355515243332, 4.141763178563418, 33.13410542850742, 10.354407946408552, 369.8956520658891, 924.739130164726, 198.15838503529784, 602.1897721340678, 2064.650647316792, 602.1897721340678, 6.902938630939034, 34.51469315469508, 6.902938630939034, 277.4217390494175, 951.160248169426, 277.4217390494175, 430.1355515243331, 2007.2992404468948, 802.9196961787545, 10.354407946408548, 33.13410542850741, 4.141763178563413, 198.15838503529775, 924.739130164726, 369.89565206588884, 286.7570343495563, 1835.24501983716, 1032.3253236584005, 14.496171124972, 28.992342249944016, 2.0708815892817136, 132.10559002353224, 845.475776150606, 475.5801240847151, 172.0542206097328, 1548.4879854875985, 1290.406654572997, 19.328228166629287, 22.08940361900494, 0.6902938630939015, 79.2633540141189, 713.3701861270716, 594.4751551058921, 86.02711030486641, 1147.0281373982186, 1577.163688922552, 24.85057907138053, 12.425289535690235, 39.63167700705945, 528.4223600941258, 726.5807451294239, 28.675703434955388, 630.8654755690198, 1892.5964267070576, 31.06322383922555, 13.210559002353111, 290.63229805176917, 871.8968941553063, 1360.8000339407185, 388.80000969734834, 14.953846526821021, 17.955886110386807, 1100.5070884823442, 314.43059670924146, 12.093484488816923, 1166.400029092047, 717.7846332874084, 44.86153958046334, 5.985295370128944, 14.963238425322388, 943.2917901277257, 580.4872554632118, 36.280453466450986, 986.9538707701863, 986.9538707701884, 89.72307916092608, 0.27205888046040594, 10.88235521841627, 12.242649620718316, 798.1699762619158, 798.1699762619177, 72.56090693290149, 822.4615589751578, 1196.3077221456845, 149.5384652682102, 0.8161766413812193, 14.691179544861999, 9.794119696574676, 665.1416468849321, 967.4787591053558, 120.9348448881692, 672.9230937069473, 1345.8461874138957, 224.30769790231605, 1.6323532827624392, 17.411768349466115, 7.617648652891407, 544.2068019967626, 1088.4136039935258, 181.4022673322544, 538.3384749655593, 1435.5692665748286, 314.03077706324376, 2.7205888046040823, 19.044121632228585, 5.713236489668569, 435.36544159741123, 1160.9745109264327, 253.9631742651572, 418.7077027509907, 1465.4769596284616, 418.7077027509907, 4.080883206906112, 19.58823939314927, 4.08088320690611, 338.61756568687554, 1185.16147990406, 338.61756568687554, 314.03077706324376, 1435.5692665748286, 538.338474965559, 5.713236489668569, 19.044121632228585, 2.7205888046040827, 253.9631742651572, 1160.9745109264331, 435.36544159741106, 224.30769790231616, 1345.8461874138957, 672.9230937069474, 7.617648652891404, 17.411768349466122, 1.6323532827624392, 181.4022673322545, 1088.4136039935258, 544.2068019967627, 149.5384652682102, 1196.3077221456845, 822.4615589751578, 9.794119696574676, 14.691179544861992,
0.8161766413812197, 120.93484488816917, 967.4787591053561, 665.1416468849321, 89.72307916092608, 986.9538707701879, 986.9538707701854, 12.24264962071831, 10.882355218416265, 0.27205888046040594, 72.56090693290149, 798.1699762619172, 798.1699762619153, 44.861539580463344, 717.7846332874083, 1166.400029092047, 14.963238425322388, 5.985295370128944, 36.280453466450986, 580.4872554632117, 943.2917901277257, 14.953846526821021, 388.80000969734834, 1360.8000339407185, 17.955886110386807, 12.093484488816923, 314.43059670924146, 1100.5070884823442, 728.9999991150366, 182.2499997787599, 6.074999992625334, 8.649772120669944, 1379.3153197874515, 344.8288299468642, 11.494294331562154, 637.8749992256597, 340.1999995870174, 18.224999977876085, 2.471363463048556, 7.41409038914568, 1206.9009048140251, 643.680482567478, 34.482882994686605, 552.8249993289004, 473.8499994247725, 36.449999955751686, 0.0950524408864824, 4.562517162551156, 6.273461098507834, 1045.9807841721463, 896.5549578618408, 68.96576598937229, 473.8499994247742, 583.1999992920278, 60.74999992625292, 0.28515732265944915, 6.27346109850785, 5.227884248756548, 896.5549578618442, 1103.4522558299584, 114.94294331562072, 400.9499995132689, 668.2499991887812, 91.12499988937927, 0.5703146453188945, 7.604195270918613, 4.277359839891719, 758.6234258830958, 1264.372376471826, 172.41441497343084, 334.1249995943908, 728.9999991150352, 127.57499984513092, 0.950524408864825, 8.554719679783448, 3.421887871913385, 632.1861882359135, 1379.315319787449, 241.3801809628031, 273.3749996681383, 765.449999070788, 170.09999979350852, 1.4257866132972423, 9.125034325102384, 2.6614683448215226, 517.2432449202936, 1448.281085776823, 321.84024128373864, 218.69999973451078, 777.5999990560397, 218.69999973451078, 1.9961012586161466, 9.315139206875294, 1.9961012586161466, 413.7945959362351, 1471.2696744399495, 413.7945959362351, 170.09999979350866, 765.4499990707884, 273.3749996681383, 2.6614683448215226, 9.125034325102385, 1.4257866132972423, 321.8402412837389, 1448.281085776824, 517.2432449202936, 127.57499984513088, 728.9999991150347, 334.12499959439066, 3.4218878719133823, 8.554719679783448, 0.9505244088648246, 241.380180962803, 1379.3153197874483, 632.1861882359132, 91.12499988937927, 668.2499991887806, 400.9499995132689, 4.277359839891721, 7.604195270918618, 0.5703146453188945, 172.41441497343084, 1264.3723764718252, 758.6234258830958, 60.74999992625286, 583.1999992920277, 473.849999424774, 5.227884248756548, 6.273461098507847, 0.28515732265944926, 114.94294331562061, 1103.452255829958, 896.554957861844, 36.44999995575165, 473.84999942477214, 552.8249993289004, 6.273461098507832, 4.5625171625511545, 0.09505244088648243, 68.96576598937223, 896.5549578618403, 1045.9807841721463, 18.22499997787609, 340.1999995870174, 637.8749992256597, 7.41409038914568, 2.471363463048556, 34.482882994686626, 643.680482567478, 1206.9009048140251, 6.074999992625334, 182.2499997787599, 728.9999991150363, 8.649772120669944, 11.494294331562154, 344.8288299468642, 1379.315319787451, 2.7824480221439876, 1719.367920948064, 382.08176021067965, 11.237698829725874, 0.6956120055359997, 2.4346420193759997, 1528.3270408427247, 719.2127251024559,
33.71309648917817, 0.023187066851200044, 1.2984757436671954, 2.1100230834591813, 1348.5238595671026, 1011.3928946753307, 67.42619297835554, 0.06956120055360031, 1.8085912143935865, 1.8085912143935932, 1179.9583771212237, 1258.6222689293086, 112.37698829725942, 0.1391224011071988, 2.2259584177151845, 1.5303464121791879, 1022.6305935050589, 1460.900847864368, 168.56548244588828, 0.23187066851199842, 2.55057735363198, 1.2752886768159908, 876.5405087186214, 1618.228631480527, 235.99167542424348, 0.3478060027679971, 2.782448022143982, 1.043418008303994, 741.6881227619083, 1730.6056197777868, 314.65556723232487, 0.486928403875196, 2.9215704232511857, 0.8347344066431954, 618.0734356349232, 1798.0318127561463, 404.5571578701323, 0.6492378718335974, 2.9679445569535905, 0.6492378718335975, 505.6964473376661, 1820.5072104155977, 505.6964473376661, 0.8347344066431954, 2.9215704232511865, 0.4869284038751958, 404.5571578701323, 1798.0318127561457, 618.0734356349229, 1.043418008303994, 2.782448022143981, 0.3478060027679972, 314.65556723232487, 1730.6056197777862, 741.6881227619083, 1.2752886768159897, 2.5505773536319776, 0.23187066851199833, 235.99167542424328, 1618.228631480527, 876.5405087186211, 1.5303464121791879, 2.2259584177151845, 0.13912240110719873, 168.56548244588834, 1460.900847864368, 1022.6305935050582, 1.8085912143935932, 1.8085912143935865, 0.06956120055360034, 112.37698829725942, 1258.6222689293081, 1179.9583771212233, 2.1100230834591813, 1.2984757436671954, 0.023187066851200012, 67.42619297835552, 1011.3928946753307,
1348.5238595671026, 2.4346420193759997, 0.6956120055359997, 33.71309648917812, 719.2127251024559, 1528.3270408427247, 2.7824480221439876, 11.237698829725897, 382.0817602106798, 1719.367920948064, 3712.5000003112104, 7425.00000062242, 3712.5000003112104, 4320.00000031018, 4320.000000310184, 720.0000000516966, 2160.000000155092, 5760.0000004135745, 2160.000000155092, 720.0000000516966, 4320.000000310184, 4320.00000031018, 4387.500000224434, 2925.0000001496232, 292.5000000149619, 2925.0000001496232, 4680.000000239403, 877.5000000448869, 1755.0000000897749, 5265.000000269316, 1755.0000000897749, 877.5000000448871,
4680.000000239403, 2925.000000149624, 292.5000000149619, 2925.0000001496232, 4387.500000224432, 4200.000000109247, 2100.0000000546274, 150.00000000390168, 3150.0000000819423, 3600.0000000936443, 450.000000011705, 2250.000000058526, 4500.0000001170565, 900.0000000234099, 1500.0000000390173, 4800.0000001248445, 1500.0000000390182, 900.0000000234104, 4500.00000011706, 2250.000000058526, 450.000000011705, 3600.0000000936443, 3150.0000000819423, 150.0000000039014, 2100.000000054626, 4200.000000109243, 3835.227272747661, 1534.0909090990676, 85.22727272772588, 3068.1818181981394, 2727.2727272872385, 255.68181818317842, 2386.3636363763303, 3579.545454564504, 511.3636363663561, 1789.7727272822494, 4090.9090909308557, 852.2727272772607, 1278.4090909158913, 4261.3636363862925, 1278.4090909158908, 852.2727272772607, 4090.9090909308557, 1789.7727272822513, 511.3636363663556, 3579.545454564504, 2386.3636363763294, 255.6818181831786, 2727.2727272872385, 3068.1818181981384, 85.22727272772585, 1534.090909099066, 3835.2272727476593, 3323.0769230648107, 1107.6923076882717, 50.34965034946678, 2769.2307692206828, 2013.9860139786767, 151.04895104840068, 2265.734265726014, 2718.88111887122, 302.0979020968014, 1812.5874125808148, 3222.377622365899, 503.4965034946719, 1409.7902097850776, 3524.475524462706, 755.2447552420059, 1057.3426573388108, 3625.174825161616, 1057.3426573388108, 755.2447552420057, 3524.475524462706, 1409.790209785076, 503.4965034946721, 3222.3776223659, 1812.587412580816, 302.0979020968014, 2718.8811188712198, 2265.734265726013, 151.04895104840077, 2013.986013978676, 2769.2307692206828, 50.34965034946678, 1107.6923076882717, 3323.0769230648107, 2677.4856438170027, 764.9958982334299, 29.422919162824087, 2294.987694700293, 1412.3001198155553, 88.26875748847281, 1941.9126647463877, 1941.912664746392, 176.53751497694444, 1618.2605539553283, 2353.8335330259324, 294.2291916282409, 1324.0313623270868, 2648.062724654176, 441.34378744236284, 1059.2250898616721, 2824.600239631133, 617.8813024193104, 823.8417365590786, 2883.446077956764, 823.8417365590786, 617.8813024193104, 2824.6002396311337, 1059.2250898616717, 441.343787442363, 2648.062724654176, 1324.0313623270868, 294.2291916282408, 2353.8335330259333, 1618.2605539553283, 176.53751497694444, 1941.9126647463909, 1941.9126647463859, 88.26875748847281, 1412.3001198155548, 2294.987694700293, 29.422919162824087, 764.9958982334299, 2677.4856438170027, 1905.8823529405947, 476.4705882351505, 15.882352941171696, 1667.6470588230275, 889.4117647056116, 47.6470588235153, 1445.2941176466113, 1238.8235294113827, 95.29411764702934, 1238.8235294113872, 1524.7058823524717, 158.82352941171587, 1048.2352941173235, 1747.0588235288722, 238.2352941175735, 873.5294117644366, 1905.8823529405904, 333.5294117646028, 714.7058823527218, 2001.176470587623, 444.7058823528055, 571.7647058821779, 2032.9411764699696, 571.7647058821779, 444.70588235280576, 2001.1764705876244, 714.7058823527218, 333.5294117646026, 1905.8823529405897, 873.5294117644362, 238.2352941175735, 1747.058823528871, 1048.2352941173235, 158.82352941171573, 1524.7058823524712, 1238.8235294113872, 95.29411764702927, 1238.8235294113822, 1445.2941176466113, 47.64705882351532, 889.4117647056116, 1667.6470588230275, 15.882352941171696, 476.4705882351505, 1905.882352940594, 1012.4998396291128, 224.99996436202437, 6.617646010647776, 899.9998574481012, 423.52934468145764, 19.85293803194365, 794.1175212777321, 595.5881409583012, 39.70587606388684, 694.8528311180206, 741.1763531925573, 66.17646010647817, 602.2057869689503, 860.2939813842133, 99.26469015971675, 516.1763888305285, 952.9410255332806, 138.9705662236034, 436.7646367027537, 1019.117485639759, 185.294088298138,
363.9705305856278, 1058.8233617036483, 238.2352563833205, 297.794070479151, 1072.0586537249435, 297.794070479151, 238.2352563833205, 1058.8233617036478, 363.9705305856276, 185.294088298138, 1019.1174856397588, 436.7646367027537, 138.9705662236033, 952.9410255332806, 516.1763888305284, 99.26469015971679, 860.2939813842133, 602.20578696895, 66.17646010647817, 741.1763531925571, 694.8528311180203, 39.70587606388683, 595.5881409583012, 794.1175212777321, 19.852938031943623, 423.52934468145764, 899.9998574481012, 6.61764601064779, 224.99996436202446, 1012.4998396291128 };

        public static int[] gamma_ind_lst0 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9,
9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21,
21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25,
25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36,
36, 36, 36, 36, 36, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
40, 40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43,
44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47,
47, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51,
51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55,
55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60,
60, 60, 60, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65, 65, 65, 65, 65,
65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69,
69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 74, 74, 74, 74, 74, 74, 74,
74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 78, 78,
78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82,
83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 87, 87, 87, 87, 87, 87, 87,
87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 92, 92,
92, 92, 92, 92, 92, 92, 92, 92, 92, 92, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 93, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 94, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 95, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 96, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 97, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98, 98,
98, 98, 98, 98, 98, 98, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 101, 101, 101, 101, 101, 101, 101,
101, 101, 101, 101, 101, 101, 101, 101, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 102, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 104, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 105, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106, 106,
106, 106, 106, 106, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 107, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 109, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 110, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 112, 112, 112, 112, 112, 112, 112, 112, 112,
112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 112, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 113, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 114, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 115, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 116, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 117, 118,
118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 118, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 119, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 120, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 121, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 122, 123, 123,
123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 123, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 124, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 125, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 126, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 128, 128, 128, 128, 128, 128, 128, 128, 128,
128, 128, 128, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 129, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 130, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 131, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 132, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 133, 134, 134, 134, 134, 134, 134, 134, 134, 134, 134, 134, 134, 134,
134, 134, 134, 134, 134, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 135, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 136, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 137, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 138, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 139, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140,
140, 141, 141, 141, 141, 141, 141, 141, 141, 141, 142, 142, 142, 142, 142, 142, 142, 142, 142, 142, 142, 143, 143, 143, 143, 143, 143, 143, 143, 143, 143, 143, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 144, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 145, 146, 146, 146, 146, 146, 146, 146, 146, 146, 146, 146, 146, 146, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 147, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 148, 149, 149, 149, 149, 149, 149, 149, 149,
149, 149, 149, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 150, 151, 151, 151, 151, 151, 151, 151, 151, 151, 151, 151, 151, 151, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 152, 153, 153, 153, 153, 153, 153, 153, 153, 153, 153, 153, 153, 153, 153, 153, 154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 154, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 155, 156, 156, 156, 156, 156, 156, 156, 156, 156, 156, 156, 157, 157, 157, 157, 157, 157, 157,
157, 157, 157, 157, 158, 158, 158, 158, 158, 158, 158, 158, 158, 158, 158, 158, 158, 159, 159, 159, 159, 159, 159, 159, 159, 159, 159, 159, 159, 159, 159, 159, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 160, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 161, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 162, 163, 163, 163, 163, 163, 163, 163, 163, 163, 163, 163, 163, 163, 163, 163, 164, 164, 164, 164, 164, 164, 164, 164, 164, 164, 164,
164, 164, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 165, 166, 166, 166, 166, 166, 166, 166, 166, 166, 166, 166, 167, 167, 167, 167, 167, 167, 167, 167, 167, 167, 167, 167, 167, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 168, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 169, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 170, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 171, 172, 172, 172, 172, 172, 172,
172, 172, 172, 172, 172, 172, 172, 172, 172, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 173, 174, 174, 174, 174, 174, 174, 174, 174, 174, 174, 174, 174, 174, 174, 174, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 175, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 176, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 177, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 178, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 179, 180,
180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 180, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 181, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 182, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 183, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 185, 186, 186, 186, 186, 186, 186, 186, 186, 186, 186, 186, 186, 186, 186,
186, 187, 187, 187, 187, 187, 187, 187, 187, 187, 187, 187, 187, 187, 187, 187, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 188, 189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 189, 190, 190, 190, 190, 190, 190, 190, 190, 191, 191, 191, 191, 191, 191, 191, 191, 191, 191, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 193, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 194, 195, 195, 195, 195, 195, 195, 195, 195, 195,
195, 195, 195, 196, 196, 196, 196, 196, 196, 196, 196, 196, 196, 196, 196, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 198, 198, 198, 198, 198, 198, 198, 198, 198, 198, 198, 198, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 199, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 201, 201, 201, 201, 201, 201, 201, 201, 201, 201, 201, 201, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 202, 203, 203, 203, 203, 203, 203, 203, 203, 203, 203, 204, 204, 204, 204, 204, 204,
204, 204, 205, 205, 205, 205, 205, 205, 206, 206, 206, 206, 206, 206, 206, 207, 207, 207, 207, 207, 207, 207, 208, 208, 208, 208, 208, 208, 208, 209, 209, 209, 209, 209, 209, 209, 210, 210, 210, 210, 210, 210, 210, 210, 211, 211, 211, 211, 211, 211, 211, 211, 211, 212, 212, 212, 212, 212, 212, 212, 212, 213, 213, 213, 213, 213, 213, 213, 214, 214, 214, 214, 214, 214, 214, 215, 215, 215, 215, 215, 215, 215, 215, 216, 216, 216, 216, 216, 216, 216, 216, 216, 217, 217, 217, 217, 217, 217, 217, 217, 217, 218, 218,
218, 218, 218, 218, 218, 218, 218, 219, 219, 219, 219, 219, 219, 219, 219, 220, 220, 220, 220, 220, 220, 220, 221, 221, 221, 221, 221, 221, 221, 222, 222, 222, 222, 222, 222, 222, 222, 223, 223, 223, 223, 223, 223, 223, 223, 223, 224, 224, 224, 224, 224, 224, 224, 224, 224, 225, 225, 225, 225, 225, 225, 225, 225, 225, 226, 226, 226, 226, 226, 226, 226, 226, 226, 227, 227, 227, 227, 227, 227, 227, 227, 227, 228, 228, 228, 228, 228, 228, 228, 228, 229, 229, 229, 229, 229, 229, 229, 230, 230, 230, 230, 230, 230,
230, 231, 231, 231, 231, 231, 231, 231, 231, 232, 232, 232, 232, 232, 232, 232, 232, 232, 233, 233, 233, 233, 233, 233, 233, 233, 233, 234, 234, 234, 234, 234, 234, 234, 234, 234, 235, 235, 235, 235, 235, 235, 235, 235, 235, 236, 236, 236, 236, 236, 236, 236, 236, 236, 237, 237, 237, 237, 237, 237, 237, 237, 237, 238, 238, 238, 238, 238, 238, 238, 238, 238, 239, 239, 239, 239, 239, 239, 239, 239, 240, 240, 240, 240, 240, 240, 240, 241, 241, 241, 241, 241, 241, 241, 242, 242, 242, 242, 242, 242, 242, 242, 243,
243, 243, 243, 243, 243, 243, 243, 243, 244, 244, 244, 244, 244, 244, 244, 244, 244, 245, 245, 245, 245, 245, 245, 245, 245, 245, 246, 246, 246, 246, 246, 246, 246, 246, 246, 247, 247, 247, 247, 247, 247, 247, 247, 247, 248, 248, 248, 248, 248, 248, 248, 248, 248, 249, 249, 249, 249, 249, 249, 249, 249, 249, 250, 250, 250, 250, 250, 250, 250, 250, 250, 251, 251, 251, 251, 251, 251, 251, 251, 251, 252, 252, 252, 252, 252, 252, 252, 252, 253, 253, 253, 253, 253, 253, 253, 254, 254, 254, 254, 254, 254, 254, 255,
255, 255, 255, 255, 255, 255, 255, 256, 256, 256, 256, 256, 256, 256, 256, 256, 257, 257, 257, 257, 257, 257, 257, 257, 257, 258, 258, 258, 258, 258, 258, 258, 258, 258, 259, 259, 259, 259, 259, 259, 259, 259, 259, 260, 260, 260, 260, 260, 260, 260, 260, 260, 261, 261, 261, 261, 261, 261, 261, 261, 261, 262, 262, 262, 262, 262, 262, 262, 262, 262, 263, 263, 263, 263, 263, 263, 263, 263, 263, 264, 264, 264, 264, 264, 264, 264, 264, 264, 265, 265, 265, 265, 265, 265, 265, 265, 265, 266, 266, 266, 266, 266, 266,
266, 266, 266, 267, 267, 267, 267, 267, 267, 267, 267, 268, 268, 268, 268, 268, 268, 268, 269, 269, 269, 269, 270, 270, 270, 270, 270, 271, 271, 271, 271, 271, 271, 272, 272, 272, 272, 272, 272, 273, 273, 273, 273, 273, 273, 274, 274, 274, 274, 274, 274, 275, 275, 275, 275, 275, 275, 276, 276, 276, 276, 276, 276, 277, 277, 277, 277, 277, 277, 278, 278, 278, 278, 278, 278, 279, 279, 279, 279, 279, 279, 280, 280, 280, 280, 280, 280, 281, 281, 281, 281, 281, 281, 282, 282, 282, 282, 282, 282, 283, 283, 283, 283,
283, 283, 284, 284, 284, 284, 284, 285, 285, 285, 285, 286, 286, 286, 287, 287, 287, 288, 288, 288, 289, 289, 289, 290, 290, 290, 291, 291, 291, 292, 292, 292, 293, 293, 293, 294, 294, 294, 295, 295, 295, 296, 296, 296, 297, 297, 297, 298, 298, 298, 299, 299, 299, 300, 300, 300, 301, 301, 301, 302, 302, 302, 303, 303, 303, 304, 304, 304, 305, 305, 305, 306, 306, 306, 307, 307, 307, 308, 308, 308, 309, 309, 309, 310, 310, 310, 311, 311, 311, 312, 312, 312, 313, 313, 313, 314, 314, 314, 315, 315, 315, 316, 316,
316, 317, 317, 317, 318, 318, 318, 319, 319, 319, 320, 320, 320, 321, 321, 321, 322, 322, 322, 323, 323, 323, 324, 324, 324, 325, 325, 325, 326, 326, 326, 327, 327, 327, 328, 328, 328, 329, 329, 329, 330, 330, 330, 331, 331, 331, 332, 332, 332, 333, 333, 333, 334, 334, 334, 335, 335, 335, 336, 336, 336, 337, 337, 337, 338, 338, 338, 339, 339, 339, 340, 340, 340, 341, 341, 341, 342, 342, 342, 343, 343, 343, 344, 344, 344, 345, 345, 345, 346, 346, 346, 347, 347, 347, 348, 348, 348, 349, 349, 349, 350, 350, 350,
351, 351, 351, 352, 352, 352, 353, 353, 353, 354, 354, 354, 355, 355, 355, 356, 356, 356, 357, 357, 357, 358, 358, 358, 359, 359, 359, 360, 360, 360, 361, 361, 361, 362, 362, 362, 363, 363, 363, 364, 364, 364, 365, 365, 365, 366, 366, 366 };

        public static int[] delta_ind_lst0 = {3, 4, 5, 7, 8, 9, 16, 17, 18, 32, 33, 34, 57, 58, 59, 93, 94, 95, 142, 143, 144, 206, 207, 208, 287, 288, 289, 3, 4, 5, 7, 8, 9, 16, 17, 18, 32, 33, 34, 57, 58, 59, 93, 94, 95, 142, 143, 144, 206, 207, 208, 287, 288, 289, 6, 10, 11, 12, 15, 19, 20, 21, 31, 35, 36, 37, 56, 60, 61, 62, 92, 96, 97, 98, 141, 145, 146, 147, 205, 209, 210, 211, 286, 290, 291, 292, 6, 11, 12, 13, 15, 20, 21, 22, 31, 36, 37, 38, 56, 61, 62, 63, 92, 97, 98, 99, 141, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 6, 12, 13, 14, 15, 21, 22, 23, 31, 37, 38, 39, 56, 62, 63, 64, 92, 98, 99, 100, 141, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 7, 8, 9, 16, 17, 18, 32, 33, 34, 57, 58, 59, 93, 94, 95, 142, 143, 144, 206, 207, 208, 287, 288, 289, 10, 11, 12, 15, 19, 20, 21, 31, 35, 36, 37, 56, 60, 61, 62, 92, 96, 97, 98, 141,
145, 146, 147, 205, 209, 210, 211, 286, 290, 291, 292, 11, 12, 13, 15, 20, 21, 22, 31, 36, 37, 38, 56,
61, 62, 63, 92, 97, 98, 99, 141, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 12, 13, 14, 15, 21, 22, 23, 31, 37, 38, 39, 56, 62, 63, 64, 92, 98, 99, 100, 141, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 16, 24, 25, 26, 32, 40, 41, 42, 57, 65, 66, 67, 93, 101, 102, 103, 142, 150, 151, 152, 206, 214, 215, 216, 287, 295, 296, 297, 16, 17, 25, 26, 27, 32, 33, 41, 42, 43, 57, 58, 66, 67, 68,
93, 94, 102, 103, 104, 142, 143, 151, 152, 153, 206, 207, 215, 216, 217, 287, 288, 296, 297, 298, 16, 17, 18, 26, 27, 28, 32, 33, 34, 42, 43, 44, 57, 58, 59, 67, 68, 69, 93, 94, 95, 103, 104, 105, 142, 143, 144, 152, 153, 154, 206, 207, 208, 216, 217, 218, 287, 288, 289, 297, 298, 299, 17, 18, 27, 28, 29, 33, 34, 43, 44, 45, 58, 59, 68, 69, 70, 94, 95, 104, 105, 106, 143, 144, 153, 154, 155, 207, 208, 217, 218, 219, 288, 289, 298, 299, 300, 18, 28, 29, 30, 34, 44, 45, 46, 59, 69, 70, 71, 95, 105, 106, 107, 144, 154, 155, 156, 208, 218, 219, 220, 289, 299, 300, 301, 16, 17, 18, 32, 33, 34, 57, 58, 59, 93, 94, 95, 142, 143, 144, 206, 207, 208, 287, 288, 289, 19, 20, 21, 31, 35, 36, 37, 56, 60, 61, 62, 92, 96, 97,
98, 141, 145, 146, 147, 205, 209, 210, 211, 286, 290, 291, 292, 20, 21, 22, 31, 36, 37, 38, 56, 61, 62, 63, 92, 97, 98, 99, 141, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 21, 22, 23, 31, 37, 38, 39, 56, 62, 63, 64, 92, 98, 99, 100, 141, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 24,
25, 26, 32, 40, 41, 42, 57, 65, 66, 67, 93, 101, 102, 103, 142, 150, 151, 152, 206, 214, 215, 216, 287, 295, 296, 297, 25, 26, 27, 32, 33, 41, 42, 43, 57, 58, 66, 67, 68, 93, 94, 102, 103, 104, 142, 143, 151, 152, 153, 206, 207, 215, 216, 217, 287, 288, 296, 297, 298, 26, 27, 28, 32, 33, 34, 42, 43, 44, 57,
58, 59, 67, 68, 69, 93, 94, 95, 103, 104, 105, 142, 143, 144, 152, 153, 154, 206, 207, 208, 216, 217, 218, 287, 288, 289, 297, 298, 299, 27, 28, 29, 33, 34, 43, 44, 45, 58, 59, 68, 69, 70, 94, 95, 104, 105, 106, 143, 144, 153, 154, 155, 207, 208, 217, 218, 219, 288, 289, 298, 299, 300, 28, 29, 30, 34, 44, 45, 46, 59, 69, 70, 71, 95, 105, 106, 107, 144, 154, 155, 156, 208, 218, 219, 220, 289, 299, 300, 301, 35, 47, 48, 49, 60, 72, 73, 74, 96, 108, 109, 110, 145, 157, 158, 159, 209, 221, 222, 223, 290, 302, 303, 304, 35, 36, 48, 49, 50, 60, 61, 73, 74, 75, 96, 97, 109, 110, 111, 145, 146, 158, 159, 160, 209, 210, 222, 223, 224, 290, 291, 303, 304, 305, 35, 36, 37, 49, 50, 51, 60, 61, 62, 74, 75, 76, 96, 97, 98, 110, 111, 112, 145, 146, 147, 159, 160, 161, 209, 210, 211, 223, 224, 225, 290, 291, 292, 304, 305, 306,
36, 37, 38, 50, 51, 52, 61, 62, 63, 75, 76, 77, 97, 98, 99, 111, 112, 113, 146, 147, 148, 160, 161, 162, 210, 211, 212, 224, 225, 226, 291, 292, 293, 305, 306, 307, 37, 38, 39, 51, 52, 53, 62, 63, 64, 76, 77, 78, 98, 99, 100, 112, 113, 114, 147, 148, 149, 161, 162, 163, 211, 212, 213, 225, 226, 227, 292, 293, 294, 306, 307, 308, 38, 39, 52, 53, 54, 63, 64, 77, 78, 79, 99, 100, 113, 114, 115, 148, 149, 162, 163, 164, 212, 213, 226, 227, 228, 293, 294, 307, 308, 309, 39, 53, 54, 55, 64, 78, 79, 80, 100, 114, 115, 116, 149, 163, 164, 165, 213, 227, 228, 229, 294, 308, 309, 310, 32, 33, 34, 57, 58, 59, 93, 94, 95,
142, 143, 144, 206, 207, 208, 287, 288, 289, 35, 36, 37, 56, 60, 61, 62, 92, 96, 97, 98, 141, 145, 146, 147, 205, 209, 210, 211, 286, 290, 291, 292, 36, 37, 38, 56, 61, 62, 63, 92, 97, 98, 99, 141, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 37, 38, 39, 56, 62, 63, 64, 92, 98, 99, 100, 141, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 40, 41, 42, 57, 65, 66, 67, 93, 101, 102, 103, 142, 150, 151, 152, 206, 214, 215, 216, 287, 295, 296, 297, 41, 42, 43, 57, 58, 66, 67, 68, 93, 94, 102, 103,
104, 142, 143, 151, 152, 153, 206, 207, 215, 216, 217, 287, 288, 296, 297, 298, 42, 43, 44, 57, 58, 59, 67, 68, 69, 93, 94, 95, 103, 104, 105, 142, 143, 144, 152, 153, 154, 206, 207, 208, 216, 217, 218, 287, 288, 289, 297, 298, 299, 43, 44, 45, 58, 59, 68, 69, 70, 94, 95, 104, 105, 106, 143, 144, 153, 154, 155, 207, 208, 217, 218, 219, 288, 289, 298, 299, 300, 44, 45, 46, 59, 69, 70, 71, 95, 105, 106, 107, 144, 154, 155, 156, 208, 218, 219, 220, 289, 299, 300, 301, 47, 48, 49, 60, 72, 73, 74, 96, 108, 109, 110, 145, 157, 158, 159, 209, 221, 222, 223, 290, 302, 303, 304, 48, 49, 50, 60, 61, 73, 74, 75, 96, 97, 109, 110, 111, 145, 146, 158, 159, 160, 209, 210, 222, 223, 224, 290, 291, 303, 304, 305, 49, 50, 51, 60, 61, 62, 74, 75, 76, 96, 97, 98, 110, 111, 112, 145, 146, 147, 159, 160, 161, 209, 210, 211, 223, 224, 225, 290, 291, 292, 304, 305, 306, 50, 51, 52, 61, 62, 63, 75, 76, 77, 97, 98, 99, 111, 112, 113, 146, 147, 148, 160, 161, 162, 210, 211, 212, 224, 225, 226, 291, 292, 293, 305, 306, 307, 51, 52, 53, 62, 63, 64, 76, 77, 78, 98, 99, 100, 112, 113, 114, 147, 148, 149, 161, 162, 163, 211, 212, 213, 225, 226, 227, 292, 293, 294, 306, 307, 308, 52, 53, 54, 63, 64, 77, 78, 79, 99, 100, 113, 114, 115, 148, 149, 162, 163, 164, 212, 213, 226, 227, 228, 293, 294, 307, 308, 309, 53, 54, 55, 64, 78, 79, 80, 100, 114, 115, 116, 149, 163, 164, 165, 213, 227, 228, 229, 294, 308, 309, 310, 65, 81, 82, 83, 101, 117, 118, 119,
150, 166, 167, 168, 214, 230, 231, 232, 295, 311, 312, 313, 65, 66, 82, 83, 84, 101, 102, 118, 119, 120, 150, 151, 167, 168, 169, 214, 215, 231, 232, 233, 295, 296, 312, 313, 314, 65, 66, 67, 83, 84, 85, 101, 102, 103, 119, 120, 121, 150, 151, 152, 168, 169, 170, 214, 215, 216, 232, 233, 234, 295, 296, 297,
313, 314, 315, 66, 67, 68, 84, 85, 86, 102, 103, 104, 120, 121, 122, 151, 152, 153, 169, 170, 171, 215, 216, 217, 233, 234, 235, 296, 297, 298, 314, 315, 316, 67, 68, 69, 85, 86, 87, 103, 104, 105, 121, 122, 123, 152, 153, 154, 170, 171, 172, 216, 217, 218, 234, 235, 236, 297, 298, 299, 315, 316, 317, 68, 69, 70, 86, 87, 88, 104, 105, 106, 122, 123, 124, 153, 154, 155, 171, 172, 173, 217, 218, 219, 235, 236,
237, 298, 299, 300, 316, 317, 318, 69, 70, 71, 87, 88, 89, 105, 106, 107, 123, 124, 125, 154, 155, 156, 172, 173, 174, 218, 219, 220, 236, 237, 238, 299, 300, 301, 317, 318, 319, 70, 71, 88, 89, 90, 106, 107, 124, 125, 126, 155, 156, 173, 174, 175, 219, 220, 237, 238, 239, 300, 301, 318, 319, 320, 71, 89, 90, 91, 107, 125, 126, 127, 156, 174, 175, 176, 220, 238, 239, 240, 301, 319, 320, 321, 57, 58, 59, 93, 94, 95, 142, 143, 144, 206, 207, 208, 287, 288, 289, 60, 61, 62, 92, 96, 97, 98, 141, 145, 146, 147, 205, 209, 210, 211, 286, 290, 291, 292, 61, 62, 63, 92, 97, 98, 99, 141, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 62, 63, 64, 92, 98, 99, 100, 141, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 65, 66, 67, 93, 101, 102, 103, 142, 150, 151, 152, 206, 214, 215, 216, 287, 295, 296, 297, 66,
67, 68, 93, 94, 102, 103, 104, 142, 143, 151, 152, 153, 206, 207, 215, 216, 217, 287, 288, 296, 297, 298, 67, 68, 69, 93, 94, 95, 103, 104, 105, 142, 143, 144, 152, 153, 154, 206, 207, 208, 216, 217, 218, 287, 288, 289, 297, 298, 299, 68, 69, 70, 94, 95, 104, 105, 106, 143, 144, 153, 154, 155, 207, 208, 217, 218, 219, 288, 289, 298, 299, 300, 69, 70, 71, 95, 105, 106, 107, 144, 154, 155, 156, 208, 218, 219, 220, 289, 299, 300, 301, 72, 73, 74, 96, 108, 109, 110, 145, 157, 158, 159, 209, 221, 222, 223, 290, 302, 303, 304, 73, 74, 75, 96, 97, 109, 110, 111, 145, 146, 158, 159, 160, 209, 210, 222, 223, 224, 290, 291, 303, 304, 305, 74, 75, 76, 96, 97, 98, 110, 111, 112, 145, 146, 147, 159, 160, 161, 209, 210, 211,
223, 224, 225, 290, 291, 292, 304, 305, 306, 75, 76, 77, 97, 98, 99, 111, 112, 113, 146, 147, 148, 160, 161, 162, 210, 211, 212, 224, 225, 226, 291, 292, 293, 305, 306, 307, 76, 77, 78, 98, 99, 100, 112, 113, 114, 147, 148, 149, 161, 162, 163, 211, 212, 213, 225, 226, 227, 292, 293, 294, 306, 307, 308, 77, 78, 79, 99, 100, 113, 114, 115, 148, 149, 162, 163, 164, 212, 213, 226, 227, 228, 293, 294, 307, 308, 309, 78, 79, 80, 100, 114, 115, 116, 149, 163, 164, 165, 213, 227, 228, 229, 294, 308, 309, 310, 81, 82,
83, 101, 117, 118, 119, 150, 166, 167, 168, 214, 230, 231, 232, 295, 311, 312, 313, 82, 83, 84, 101, 102, 118, 119, 120, 150, 151, 167, 168, 169, 214, 215, 231, 232, 233, 295, 296, 312, 313, 314, 83, 84, 85, 101, 102, 103, 119, 120, 121, 150, 151, 152, 168, 169, 170, 214, 215, 216, 232, 233, 234, 295, 296, 297, 313, 314, 315, 84, 85, 86, 102, 103, 104, 120, 121, 122, 151, 152, 153, 169, 170, 171, 215, 216, 217, 233, 234, 235, 296, 297, 298, 314, 315, 316, 85, 86, 87, 103, 104, 105, 121, 122, 123, 152, 153, 154, 170, 171, 172, 216, 217, 218, 234, 235, 236, 297, 298, 299, 315, 316, 317, 86, 87, 88, 104, 105, 106, 122, 123, 124, 153, 154, 155, 171, 172, 173, 217, 218, 219, 235, 236, 237, 298, 299, 300, 316, 317, 318, 87, 88, 89, 105, 106, 107, 123, 124, 125, 154, 155, 156, 172, 173, 174, 218, 219, 220, 236, 237, 238, 299, 300, 301, 317, 318, 319, 88, 89, 90, 106, 107, 124, 125, 126, 155, 156, 173, 174, 175, 219, 220, 237, 238, 239, 300, 301, 318, 319, 320, 89, 90, 91, 107, 125, 126, 127, 156, 174, 175, 176, 220, 238,
239, 240, 301, 319, 320, 321, 108, 128, 129, 130, 157, 177, 178, 179, 221, 241, 242, 243, 302, 322, 323, 324, 108, 109, 129, 130, 131, 157, 158, 178, 179, 180, 221, 222, 242, 243, 244, 302, 303, 323, 324, 325, 108, 109, 110, 130, 131, 132, 157, 158, 159, 179, 180, 181, 221, 222, 223, 243, 244, 245, 302, 303, 304, 324, 325, 326, 109, 110, 111, 131, 132, 133, 158, 159, 160, 180, 181, 182, 222, 223, 224, 244, 245, 246, 303, 304, 305, 325, 326, 327, 110, 111, 112, 132, 133, 134, 159, 160, 161, 181, 182, 183, 223,
224, 225, 245, 246, 247, 304, 305, 306, 326, 327, 328, 111, 112, 113, 133, 134, 135, 160, 161, 162, 182, 183, 184, 224, 225, 226, 246, 247, 248, 305, 306, 307, 327, 328, 329, 112, 113, 114, 134, 135, 136, 161, 162, 163, 183, 184, 185, 225, 226, 227, 247, 248, 249, 306, 307, 308, 328, 329, 330, 113, 114, 115, 135, 136, 137, 162, 163, 164, 184, 185, 186, 226, 227, 228, 248, 249, 250, 307, 308, 309, 329, 330, 331, 114, 115, 116, 136, 137, 138, 163, 164, 165, 185, 186, 187, 227, 228, 229, 249, 250, 251, 308, 309,
310, 330, 331, 332, 115, 116, 137, 138, 139, 164, 165, 186, 187, 188, 228, 229, 250, 251, 252, 309, 310, 331, 332, 333, 116, 138, 139, 140, 165, 187, 188, 189, 229, 251, 252, 253, 310, 332, 333, 334, 93, 94, 95, 142, 143, 144, 206, 207, 208, 287, 288, 289, 96, 97, 98, 141, 145, 146, 147, 205, 209, 210, 211,
286, 290, 291, 292, 97, 98, 99, 141, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 98, 99, 100, 141, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 101, 102, 103, 142, 150, 151, 152, 206, 214, 215, 216, 287, 295, 296, 297, 102, 103, 104, 142, 143, 151, 152, 153, 206, 207, 215, 216, 217, 287, 288, 296, 297, 298, 103, 104, 105, 142, 143, 144, 152, 153, 154, 206, 207, 208, 216, 217, 218, 287, 288, 289, 297, 298, 299, 104, 105, 106, 143, 144, 153, 154, 155, 207, 208, 217, 218, 219, 288, 289, 298,
299, 300, 105, 106, 107, 144, 154, 155, 156, 208, 218, 219, 220, 289, 299, 300, 301, 108, 109, 110, 145, 157, 158, 159, 209, 221, 222, 223, 290, 302, 303, 304, 109, 110, 111, 145, 146, 158, 159, 160, 209, 210, 222, 223, 224, 290, 291, 303, 304, 305, 110, 111, 112, 145, 146, 147, 159, 160, 161, 209, 210, 211, 223, 224, 225, 290, 291, 292, 304, 305, 306, 111, 112, 113, 146, 147, 148, 160, 161, 162, 210, 211, 212, 224, 225, 226, 291, 292, 293, 305, 306, 307, 112, 113, 114, 147, 148, 149, 161, 162, 163, 211, 212,
213, 225, 226, 227, 292, 293, 294, 306, 307, 308, 113, 114, 115, 148, 149, 162, 163, 164, 212, 213, 226, 227, 228, 293, 294, 307, 308, 309, 114, 115, 116, 149, 163, 164, 165, 213, 227, 228, 229, 294, 308, 309, 310, 117, 118, 119, 150, 166, 167, 168, 214, 230, 231, 232, 295, 311, 312, 313, 118, 119, 120, 150, 151, 167, 168, 169, 214, 215, 231, 232, 233, 295, 296, 312, 313, 314, 119, 120, 121, 150, 151, 152, 168, 169, 170, 214, 215, 216, 232, 233, 234, 295, 296, 297, 313, 314, 315, 120, 121, 122, 151, 152, 153,
169, 170, 171, 215, 216, 217, 233, 234, 235, 296, 297, 298, 314, 315, 316, 121, 122, 123, 152, 153, 154, 170, 171, 172, 216, 217, 218, 234, 235, 236, 297, 298, 299, 315, 316, 317, 122, 123, 124, 153, 154, 155, 171, 172, 173, 217, 218, 219, 235, 236, 237, 298, 299, 300, 316, 317, 318, 123, 124, 125, 154, 155, 156, 172, 173, 174, 218, 219, 220, 236, 237, 238, 299, 300, 301, 317, 318, 319, 124, 125, 126, 155, 156, 173, 174, 175, 219, 220, 237, 238, 239, 300, 301, 318, 319, 320, 125, 126, 127, 156, 174, 175, 176,
220, 238, 239, 240, 301, 319, 320, 321, 128, 129, 130, 157, 177, 178, 179, 221, 241, 242, 243, 302, 322, 323, 324, 129, 130, 131, 157, 158, 178, 179, 180, 221, 222, 242, 243, 244, 302, 303, 323, 324, 325, 130, 131, 132, 157, 158, 159, 179, 180, 181, 221, 222, 223, 243, 244, 245, 302, 303, 304, 324, 325, 326, 131, 132, 133, 158, 159, 160, 180, 181, 182, 222, 223, 224, 244, 245, 246, 303, 304, 305, 325, 326, 327, 132, 133, 134, 159, 160, 161, 181, 182, 183, 223, 224, 225, 245, 246, 247, 304, 305, 306, 326, 327,
328, 133, 134, 135, 160, 161, 162, 182, 183, 184, 224, 225, 226, 246, 247, 248, 305, 306, 307, 327, 328, 329, 134, 135, 136, 161, 162, 163, 183, 184, 185, 225, 226, 227, 247, 248, 249, 306, 307, 308, 328, 329, 330, 135, 136, 137, 162, 163, 164, 184, 185, 186, 226, 227, 228, 248, 249, 250, 307, 308, 309, 329, 330, 331, 136, 137, 138, 163, 164, 165, 185, 186, 187, 227, 228, 229, 249, 250, 251, 308, 309, 310, 330, 331, 332, 137, 138, 139, 164, 165, 186, 187, 188, 228, 229, 250, 251, 252, 309, 310, 331, 332, 333,
138, 139, 140, 165, 187, 188, 189, 229, 251, 252, 253, 310, 332, 333, 334, 166, 190, 191, 192, 230, 254, 255, 256, 311, 335, 336, 337, 166, 167, 191, 192, 193, 230, 231, 255, 256, 257, 311, 312, 336, 337, 338, 166, 167, 168, 192, 193, 194, 230, 231, 232, 256, 257, 258, 311, 312, 313, 337, 338, 339, 167, 168, 169, 193, 194, 195, 231, 232, 233, 257, 258, 259, 312, 313, 314, 338, 339, 340, 168, 169, 170, 194, 195, 196, 232, 233, 234, 258, 259, 260, 313, 314, 315, 339, 340, 341, 169, 170, 171, 195, 196, 197, 233,
234, 235, 259, 260, 261, 314, 315, 316, 340, 341, 342, 170, 171, 172, 196, 197, 198, 234, 235, 236, 260, 261, 262, 315, 316, 317, 341, 342, 343, 171, 172, 173, 197, 198, 199, 235, 236, 237, 261, 262, 263, 316, 317, 318, 342, 343, 344, 172, 173, 174, 198, 199, 200, 236, 237, 238, 262, 263, 264, 317, 318, 319, 343, 344, 345, 173, 174, 175, 199, 200, 201, 237, 238, 239, 263, 264, 265, 318, 319, 320, 344, 345, 346, 174, 175, 176, 200, 201, 202, 238, 239, 240, 264, 265, 266, 319, 320, 321, 345, 346, 347, 175, 176,
201, 202, 203, 239, 240, 265, 266, 267, 320, 321, 346, 347, 348, 176, 202, 203, 204, 240, 266, 267, 268, 321, 347, 348, 349, 142, 143, 144, 206, 207, 208, 287, 288, 289, 145, 146, 147, 205, 209, 210, 211, 286, 290, 291, 292, 146, 147, 148, 205, 210, 211, 212, 286, 291, 292, 293, 147, 148, 149, 205, 211, 212, 213, 286, 292, 293, 294, 150, 151, 152, 206, 214, 215, 216, 287, 295, 296, 297, 151, 152, 153, 206, 207, 215, 216, 217, 287, 288, 296, 297, 298, 152, 153, 154, 206, 207, 208, 216, 217, 218, 287, 288, 289,
297, 298, 299, 153, 154, 155, 207, 208, 217, 218, 219, 288, 289, 298, 299, 300, 154, 155, 156, 208, 218, 219, 220, 289, 299, 300, 301, 157, 158, 159, 209, 221, 222, 223, 290, 302, 303, 304, 158, 159, 160, 209, 210, 222, 223, 224, 290, 291, 303, 304, 305, 159, 160, 161, 209, 210, 211, 223, 224, 225, 290, 291, 292, 304, 305, 306, 160, 161, 162, 210, 211, 212, 224, 225, 226, 291, 292, 293, 305, 306, 307, 161, 162, 163, 211, 212, 213, 225, 226, 227, 292, 293, 294, 306, 307, 308, 162, 163, 164, 212, 213, 226, 227,
228, 293, 294, 307, 308, 309, 163, 164, 165, 213, 227, 228, 229, 294, 308, 309, 310, 166, 167, 168, 214, 230, 231, 232, 295, 311, 312, 313, 167, 168, 169, 214, 215, 231, 232, 233, 295, 296, 312, 313, 314, 168, 169, 170, 214, 215, 216, 232, 233, 234, 295, 296, 297, 313, 314, 315, 169, 170, 171, 215, 216, 217, 233, 234, 235, 296, 297, 298, 314, 315, 316, 170, 171, 172, 216, 217, 218, 234, 235, 236, 297, 298, 299, 315, 316, 317, 171, 172, 173, 217, 218, 219, 235, 236, 237, 298, 299, 300, 316, 317, 318, 172, 173,
174, 218, 219, 220, 236, 237, 238, 299, 300, 301, 317, 318, 319, 173, 174, 175, 219, 220, 237, 238, 239, 300, 301, 318, 319, 320, 174, 175, 176, 220, 238, 239, 240, 301, 319, 320, 321, 177, 178, 179, 221, 241, 242, 243, 302, 322, 323, 324, 178, 179, 180, 221, 222, 242, 243, 244, 302, 303, 323, 324, 325, 179, 180, 181, 221, 222, 223, 243, 244, 245, 302, 303, 304, 324, 325, 326, 180, 181, 182, 222, 223, 224, 244, 245, 246, 303, 304, 305, 325, 326, 327, 181, 182, 183, 223, 224, 225, 245, 246, 247, 304, 305, 306,
326, 327, 328, 182, 183, 184, 224, 225, 226, 246, 247, 248, 305, 306, 307, 327, 328, 329, 183, 184, 185, 225, 226, 227, 247, 248, 249, 306, 307, 308, 328, 329, 330, 184, 185, 186, 226, 227, 228, 248, 249, 250, 307, 308, 309, 329, 330, 331, 185, 186, 187, 227, 228, 229, 249, 250, 251, 308, 309, 310, 330, 331, 332, 186, 187, 188, 228, 229, 250, 251, 252, 309, 310, 331, 332, 333, 187, 188, 189, 229, 251, 252, 253, 310, 332, 333, 334, 190, 191, 192, 230, 254, 255, 256, 311, 335, 336, 337, 191, 192, 193, 230, 231,
255, 256, 257, 311, 312, 336, 337, 338, 192, 193, 194, 230, 231, 232, 256, 257, 258, 311, 312, 313, 337, 338, 339, 193, 194, 195, 231, 232, 233, 257, 258, 259, 312, 313, 314, 338, 339, 340, 194, 195, 196, 232, 233, 234, 258, 259, 260, 313, 314, 315, 339, 340, 341, 195, 196, 197, 233, 234, 235, 259, 260, 261, 314, 315, 316, 340, 341, 342, 196, 197, 198, 234, 235, 236, 260, 261, 262, 315, 316, 317, 341, 342, 343, 197, 198, 199, 235, 236, 237, 261, 262, 263, 316, 317, 318, 342, 343, 344, 198, 199, 200, 236, 237,
238, 262, 263, 264, 317, 318, 319, 343, 344, 345, 199, 200, 201, 237, 238, 239, 263, 264, 265, 318, 319, 320, 344, 345, 346, 200, 201, 202, 238, 239, 240, 264, 265, 266, 319, 320, 321, 345, 346, 347, 201, 202, 203, 239, 240, 265, 266, 267, 320, 321, 346, 347, 348, 202, 203, 204, 240, 266, 267, 268, 321, 347, 348, 349, 241, 269, 270, 271, 322, 350, 351, 352, 241, 242, 270, 271, 272, 322, 323, 351, 352, 353, 241, 242, 243, 271, 272, 273, 322, 323, 324, 352, 353, 354, 242, 243, 244, 272, 273, 274, 323, 324, 325,
353, 354, 355, 243, 244, 245, 273, 274, 275, 324, 325, 326, 354, 355, 356, 244, 245, 246, 274, 275, 276, 325, 326, 327, 355, 356, 357, 245, 246, 247, 275, 276, 277, 326, 327, 328, 356, 357, 358, 246, 247, 248, 276, 277, 278, 327, 328, 329, 357, 358, 359, 247, 248, 249, 277, 278, 279, 328, 329, 330, 358, 359, 360, 248, 249, 250, 278, 279, 280, 329, 330, 331, 359, 360, 361, 249, 250, 251, 279, 280, 281, 330, 331, 332, 360, 361, 362, 250, 251, 252, 280, 281, 282, 331, 332, 333, 361, 362, 363, 251, 252, 253, 281,
282, 283, 332, 333, 334, 362, 363, 364, 252, 253, 282, 283, 284, 333, 334, 363, 364, 365, 253, 283, 284, 285, 334, 364, 365, 366, 206, 207, 208, 287, 288, 289, 209, 210, 211, 286, 290, 291, 292, 210, 211, 212, 286, 291, 292, 293, 211, 212, 213, 286, 292, 293, 294, 214, 215, 216, 287, 295, 296, 297, 215, 216, 217, 287, 288, 296, 297, 298, 216, 217, 218, 287, 288, 289, 297, 298, 299, 217, 218, 219, 288, 289, 298, 299, 300, 218, 219, 220, 289, 299, 300, 301, 221, 222, 223, 290, 302, 303, 304, 222, 223, 224, 290,
291, 303, 304, 305, 223, 224, 225, 290, 291, 292, 304, 305, 306, 224, 225, 226, 291, 292, 293, 305, 306, 307, 225, 226, 227, 292, 293, 294, 306, 307, 308, 226, 227, 228, 293, 294, 307, 308, 309, 227, 228, 229, 294, 308, 309, 310, 230, 231, 232, 295, 311, 312, 313, 231, 232, 233, 295, 296, 312, 313, 314, 232, 233, 234, 295, 296, 297, 313, 314, 315, 233, 234, 235, 296, 297, 298, 314, 315, 316, 234, 235, 236, 297, 298, 299, 315, 316, 317, 235, 236, 237, 298, 299, 300, 316, 317, 318, 236, 237, 238, 299, 300, 301,
317, 318, 319, 237, 238, 239, 300, 301, 318, 319, 320, 238, 239, 240, 301, 319, 320, 321, 241, 242, 243, 302, 322, 323, 324, 242, 243, 244, 302, 303, 323, 324, 325, 243, 244, 245, 302, 303, 304, 324, 325, 326, 244, 245, 246, 303, 304, 305, 325, 326, 327, 245, 246, 247, 304, 305, 306, 326, 327, 328, 246, 247, 248, 305, 306, 307, 327, 328, 329, 247, 248, 249, 306, 307, 308, 328, 329, 330, 248, 249, 250, 307, 308, 309, 329, 330, 331, 249, 250, 251, 308, 309, 310, 330, 331, 332, 250, 251, 252, 309, 310, 331, 332,
333, 251, 252, 253, 310, 332, 333, 334, 254, 255, 256, 311, 335, 336, 337, 255, 256, 257, 311, 312, 336, 337, 338, 256, 257, 258, 311, 312, 313, 337, 338, 339, 257, 258, 259, 312, 313, 314, 338, 339, 340, 258, 259, 260, 313, 314, 315, 339, 340, 341, 259, 260, 261, 314, 315, 316, 340, 341, 342, 260, 261, 262, 315, 316, 317, 341, 342, 343, 261, 262, 263, 316, 317, 318, 342, 343, 344, 262, 263, 264, 317, 318, 319, 343, 344, 345, 263, 264, 265, 318, 319, 320, 344, 345, 346, 264, 265, 266, 319, 320, 321, 345, 346,
347, 265, 266, 267, 320, 321, 346, 347, 348, 266, 267, 268, 321, 347, 348, 349, 269, 270, 271, 322, 350, 351, 352, 270, 271, 272, 322, 323, 351, 352, 353, 271, 272, 273, 322, 323, 324, 352, 353, 354, 272, 273, 274, 323, 324, 325, 353, 354, 355, 273, 274, 275, 324, 325, 326, 354, 355, 356, 274, 275, 276, 325, 326, 327, 355, 356, 357, 275, 276, 277, 326, 327, 328, 356, 357, 358, 276, 277, 278, 327, 328, 329, 357, 358, 359, 277, 278, 279, 328, 329, 330, 358, 359, 360, 278, 279, 280, 329, 330, 331, 359, 360, 361,
279, 280, 281, 330, 331, 332, 360, 361, 362, 280, 281, 282, 331, 332, 333, 361, 362, 363, 281, 282, 283, 332, 333, 334, 362, 363, 364, 282, 283, 284, 333, 334, 363, 364, 365, 283, 284, 285, 334, 364, 365, 366, 335, 367, 368, 369, 335, 336, 368, 369, 370, 335, 336, 337, 369, 370, 371, 336, 337, 338, 370, 371, 372, 337, 338, 339, 371, 372, 373, 338, 339, 340, 372, 373, 374, 339, 340, 341, 373, 374, 375, 340, 341, 342, 374, 375, 376, 341, 342, 343, 375, 376, 377, 342, 343, 344, 376, 377, 378, 343, 344, 345, 377,
378, 379, 344, 345, 346, 378, 379, 380, 345, 346, 347, 379, 380, 381, 346, 347, 348, 380, 381, 382, 347, 348, 349, 381, 382, 383, 348, 349, 382, 383, 384, 349, 383, 384, 385, 287, 288, 289, 290, 291, 292, 291, 292, 293, 292, 293, 294, 295, 296, 297, 296, 297, 298, 297, 298, 299, 298, 299, 300, 299, 300, 301, 302, 303, 304, 303, 304, 305, 304, 305, 306, 305, 306, 307, 306, 307, 308, 307, 308, 309, 308, 309, 310, 311, 312, 313, 312, 313, 314, 313, 314, 315, 314, 315, 316, 315, 316, 317, 316, 317, 318, 317, 318,
319, 318, 319, 320, 319, 320, 321, 322, 323, 324, 323, 324, 325, 324, 325, 326, 325, 326, 327, 326, 327, 328, 327, 328, 329, 328, 329, 330, 329, 330, 331, 330, 331, 332, 331, 332, 333, 332, 333, 334, 335, 336, 337, 336, 337, 338, 337, 338, 339, 338, 339, 340, 339, 340, 341, 340, 341, 342, 341, 342, 343, 342, 343, 344, 343, 344, 345, 344, 345, 346, 345, 346, 347, 346, 347, 348, 347, 348, 349, 350, 351, 352, 351, 352, 353, 352, 353, 354, 353, 354, 355, 354, 355, 356, 355, 356, 357, 356, 357, 358, 357, 358, 359,
358, 359, 360, 359, 360, 361, 360, 361, 362, 361, 362, 363, 362, 363, 364, 363, 364, 365, 364, 365, 366, 367, 368, 369, 368, 369, 370, 369, 370, 371, 370, 371, 372, 371, 372, 373, 372, 373, 374, 373, 374, 375, 374, 375, 376, 375, 376, 377, 376, 377, 378, 377, 378, 379, 378, 379, 380, 379, 380, 381, 380, 381, 382, 381, 382, 383, 382, 383, 384, 383, 384, 385};
    } //Variables imported from the hydrogen parameters. These correspond to the gamma/delta lists/indices. Call by variables.{variable name here}. {Line 361-363}
    public static class TransposeRowsColumnsExtension //Transpose method for NxM matrices. Works for any variable type. Call by {variable name}.TransposeRowsAndColumns(). {Line 618, 710, 780}
    {
        /// <summary>
        /// Transposes the rows and columns of a two-dimensional array
        /// </summary>
        /// <typeparam name="T">The type of the items in the array</typeparam>
        /// <param name="arr">The array</param>
        /// <returns>A new array with rows and columns transposed</returns>
        public static T[,] TransposeRowsAndColumns<T>(this T[,] arr)
        {
            int rowCount = arr.GetLength(0);
            int columnCount = arr.GetLength(1);
            T[,] transposed = new T[columnCount, rowCount];
            for (int column = 0; column < columnCount; column++)
            {
                for (int row = 0; row < rowCount; row++)
                {
                    transposed[column, row] = arr[row, column];
                }
            }
            return transposed;
        }

        /// <summary>
        /// Transposes the rows and columns of a jagged array
        /// </summary>
        /// <typeparam name="T">The type of the items in the array</typeparam>
        /// <param name="arr">The array</param>
        /// <returns>A new array with rows and columns transposed</returns>
        public static T[][] TransposeRowsAndColumns<T>(this T[][] arr)
        {
            int rowCount = arr.Length;
            int columnCount = arr[0].Length;
            T[][] transposed = new T[columnCount][];
            for (int column = 0; column < columnCount; column++)
            {
                transposed[column] = new T[rowCount];
                for (int row = 0; row < rowCount; row++)
                {
                    transposed[column][row] = arr[row][column];
                }
            }
            return transposed;
        }

        /// <summary>
        /// Transposes the rows and columns of a string
        /// </summary>
        /// <param name="str">The string</param>
        /// <param name="rowDelimiter">The delimiter of the rows</param>
        /// <param name="columnDelimiter">The delimiter of the columns</param>
        /// <returns>A new string with rows and columns transposed</returns>
        public static string TransposeRowsAndColumns(this string str, string rowDelimiter, string columnDelimiter)
        {
            string[] rows = str.Split(new string[] { rowDelimiter }, StringSplitOptions.None);
            string[][] arr = new string[rows.Length][];
            for (int i = 0; i < rows.Length; i++)
            {
                arr[i] = rows[i].Split(new string[] { columnDelimiter }, StringSplitOptions.None);
            }
            string[][] transposed = TransposeRowsAndColumns(arr);
            string[] transposedRows = new string[transposed.Length];
            for (int i = 0; i < transposed.Length; i++)
            {
                transposedRows[i] = String.Join(columnDelimiter, transposed[i]);
            }
            return String.Join(rowDelimiter, transposedRows);
        }
    }
    public static class matrix_mult //Performs Matrix multiplication for NxM matrices. Does NOT work for complex variable arrays. Call by matrix_mult.Multiply(variable 1, variable 2). Not referenced as of yet.
    {
        public static double[,] Multiply(double[,] matrix1, double[,] matrix2)
        {
            // cahing matrix lengths for better performance  
            var matrix1Rows = matrix1.GetLength(0);
            var matrix1Cols = matrix1.GetLength(1);
            var matrix2Rows = matrix2.GetLength(0);
            var matrix2Cols = matrix2.GetLength(1);

            // checking if product is defined  
            if (matrix1Cols != matrix2Rows)
                throw new InvalidOperationException
                  ("Product is undefined. n columns of first matrix must equal to n rows of second matrix");

            // creating the final product matrix  
            double[,] product = new double[matrix1Rows, matrix2Cols];

            // looping through matrix 1 rows  
            for (int matrix1_row = 0; matrix1_row < matrix1Rows; matrix1_row++)
            {
                // for each matrix 1 row, loop through matrix 2 columns  
                for (int matrix2_col = 0; matrix2_col < matrix2Cols; matrix2_col++)
                {
                    // loop through matrix 1 columns to calculate the dot product  
                    for (int matrix1_col = 0; matrix1_col < matrix1Cols; matrix1_col++)
                    {
                        product[matrix1_row, matrix2_col] +=
                          matrix1[matrix1_row, matrix1_col] *
                          matrix2[matrix1_col, matrix2_col];
                    }
                }
            }

            return product;
        }

    }
    public static class matrix_mult_C //Performs Matrix multiplication for NxM matrices. Use specifically for complex variable arrays. Call by matrix_mult_C.Multiply(variable 1, variable 2). {Line 486, 488, 543, 560-561, 583-584, 617-618}
    {
        public static Complex[,] Multiply(Complex[,] matrix1, Complex[,] matrix2)
        {
            // cahing matrix lengths for better performance  
            var matrix1Rows = matrix1.GetLength(0);
            var matrix1Cols = matrix1.GetLength(1);
            var matrix2Rows = matrix2.GetLength(0);
            var matrix2Cols = matrix2.GetLength(1);

            // checking if product is defined  
            if (matrix1Cols != matrix2Rows)
                throw new InvalidOperationException
                  ("Product is undefined. n columns of first matrix must equal to n rows of second matrix");

            // creating the final product matrix  
            Complex[,] product = new Complex[matrix1Rows, matrix2Cols];

            // looping through matrix 1 rows  
            for (int matrix1_row = 0; matrix1_row < matrix1Rows; matrix1_row++)
            {
                // for each matrix 1 row, loop through matrix 2 columns  
                for (int matrix2_col = 0; matrix2_col < matrix2Cols; matrix2_col++)
                {
                    // loop through matrix 1 columns to calculate the dot product  
                    for (int matrix1_col = 0; matrix1_col < matrix1Cols; matrix1_col++)
                    {
                        product[matrix1_row, matrix2_col] = Complex.Add(product[matrix1_row, matrix2_col], Complex.Multiply(matrix1[matrix1_row, matrix1_col],
                          matrix2[matrix1_col, matrix2_col]));
                    }
                }
            }

            return product;
        }

    }

}


/* References:
[1]: Abramowitz, Milton; Stegun, Irene Ann, eds. (1983) [June 1964]. "Chapter 22". 
Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical Tables. Applied Mathematics Series. 55 
(Ninth reprint with additional corrections of tenth original printing with corrections (December 1972); first ed.). 
Washington D.C.; New York: United States Department of Commerce, National Bureau of Standards; Dover Publications. p. 773. ISBN 978-0-486-61272-0. LCCN 64-60036. MR 0167642. LCCN 65-12253.
[2]: “§14.10 Recurrence Relations and Derivatives.” DLMF, dlmf.nist.gov/14.10. 
[3]: Griffiths, D. J. (2004). Introduction to Quantum Mechanics (2nd Edition). Pearson Prentice Hall. 152. ISBN: 0131118927

 */