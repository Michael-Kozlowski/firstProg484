//============================================================================
// SimplePend.cs  Defines a class for simulating a simple pendulum 
//============================================================================
using System;

namespace Sim 
{
    public class SimplePend
    {
        private double len = 1.1;     // pendulum length
        private double g = 9.81;      // gravitation field stregth 

        // int n = 2;                 // number of states
        public double[] theta;       // array of states
        public double[] w;           // right side of the equation value
        public double[] k1;           // RK slopes
        public double[] k2; 
        public double theta_answer;
        public double w_answer;
        string s;
        //--------------------------------------------------------------------
        // constructor
        //--------------------------------------------------------------------
        public SimplePend()
        {
            theta = new double[4];       // holds values for theta and w inbetween each step 
            
            w = new double[4];           // dtheta/dt and dw/dt 
 
            k1 = new double[4];        // RK slopes 

            k2 = new double[4]; 

            theta[0] = 1;

            w[0] = 0;

            s = " ";
            theta_answer = 0.0;
            w_answer = 0.0;
            // x[0] = 1.0;  // initial condition for theta in rad
            // x[1] = 0.0;  // initial condition for w (omega) in rad/s
        }

        

        //--------------------------------------------------------------------
        // RungeKutta: calculates answer
        //-------------------------------------------------------------------- 
        /*public void RungeKutta(double dt)
        {
            theta[0] = 1.0;
            w[0] = 0.0;
            RKslope(theta, w, k1, k2, dt);

            theta_answer = theta[0] + 1/6*(k1[0] + 2*k1[1] + 2*k1[2] + k1[3]) * dt;
            w_answer     = w[0]     + 1/6*(k2[0] + 2*k2[1] + 2*k2[2] + k2[3]) * dt;
        }
        */
        //--------------------------------------------------------------------
        // RKslope: Calculates the slopes (k) for each iteration 
        //--------------------------------------------------------------------

        public void RungeKutta(double dt, double t, double tEnd)   //, double tt, double ww
        {
            while (t < tEnd)
            {

            
                for (int n=0;n<4;++n)
                {

                    if (n == 0)
                    {
                    // RK4internal(theta[n], w[n], null, null);
                        k1[0] = w[0];
                        k2[0] = -g/len*Math.Sin(theta[0]);
                        theta[1] = theta[0] + k1[0]*dt/2;
                        w[1]     = w[0]     + k2[0]*dt/2;
                    }
                    else if (n == 1)
                    {
                    // RK4internal(theta[n], w[n], k1[n-1], k2[n-1]);
                        k1[1] = w[1];
                        k2[1] = -g/len*Math.Sin(theta[1]);
                        theta[2] = theta[0] + k1[1]*dt/2;
                        w[2]     = w[0]     + k2[1]*dt/2;
                    }

                    else if (n == 2)
                    {
                    // RK4internal(theta[n], w[n], k1[n-1], k2[n-1]);
                        k1[2] = w[2];
                        k2[2] = -g/len*Math.Sin(theta[2]);
                        theta[3] = theta[0] + k1[2]*dt/2;
                        w[3]     = w[0]     + k2[2]*dt/2;
                    }

                    else
                    {
                    // RK4internal(theta[n], w[n], null, null);
                        k1[3] = w[3];
                        k2[3] = -g/len*Math.Sin(theta[3]);
                    }

                }
                // Current error: Value of theta answer is being passed as 1 rather than .99875
                theta_answer = theta[0] + 1.0/6.0*(k1[0] + 2*k1[1] + 2*k1[2] + k1[3]) * dt;
                w_answer     = w[0]     + 1.0/6.0*(k2[0] + 2*k2[1] + 2*k2[2] + k2[3]) * dt;  // need to pass theta_ answer and w_answer to theta[0] and w[0]

                theta[0] = theta_answer;
                w[0] = w_answer;
                
                s = $"{t.ToString("0.####")}   {theta_answer.ToString("0.####")}   {w_answer.ToString("0.####")}";
                Console.WriteLine(s);

                t += 0.02;
            }
        }
        
        //--------------------------------------------------------------------
        //  step: perform one integration step via Euler's Method
        //--------------------------------------------------------------------
     /* public void step(double dt)
        {
            rhsFunc(x, f);
            int i; 
            for(i=0;i<n;++i)
            {
                x[i] = x[i] + f[i] * dt;
            }
            //Console.WriteLine($"{f[0].ToString()}   {f[1].ToString()}");
        }
        
        //--------------------------------------------------------------------
        // rhsFunc: Calculate the right handside of the pendulum equations 
        //--------------------------------------------------------------------
        public void rhsFunc(double[] st, double[] ff)
        {
            ff[0] = st[1];
            ff[1] = -g/len * Math.Sin(st[0]);
        }
        */
       

        //--------------------------------------------------------------------
        // Getters and Setters
        //--------------------------------------------------------------------
        public double L
        {
            get {return(len);}

            set 
            {
                if(value > 0.0)
                    len = value;
            }
        }

        public double G
        {
            get {return(g);}

            set
            {
                if (value >= 0.0)
                    g = value;
            }
        }

        public double thetaans
        {
            get{return theta_answer;}

            set {theta_answer = value;}
        }

        public double thetaDotans
        {
            get{return w_answer;}

            set {w_answer = value;}
        }
    }
}