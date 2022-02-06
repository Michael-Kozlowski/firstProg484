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

        // int n = 2;                 // number of states (Euler's Method)
        public double[] theta;        // array of states (RK Method)
        public double[] w;            // right side of the equation value (RK Method)
        public double[] k1;           // RK slopes (RK Method)
        public double[] k2;           // RK slopes (RK Method)
        public double theta_answer;   // Answer for theta value (RK Method)
        public double w_answer;       // Answer for w value (RK Method)

        //--------------------------------------------------------------------
        // constructor
        //--------------------------------------------------------------------
        public SimplePend()
        {
            theta = new double[4];       // Holds values for theta inbetween each step (RK Method)
            
            w = new double[4];           // Holds values for w inbetween each step (RK Method)
 
            k1 = new double[4];          // Holds slope values (RK Method)

            k2 = new double[4];          // Holds slope values (RK Method)

            theta[0] = 1;                // Initial condition for theta in rad (RK Method)

            w[0] = 0;                    // Initial condtion for w in rad/s (RK Method)

            // x[0] = 1.0;  // Initial condition for theta in rad (Euler's Method)
            // x[1] = 0.0;  // Initial condition for w (omega) in rad/s (Euler's Method)
        }

        //--------------------------------------------------------------------
        // RKslope: Calculates the slopes (k) for each iteration 
        //--------------------------------------------------------------------
        public void RungeKutta(double dt, double t, double tEnd)   
        {            
            for (int n=0;n<4;++n)
            {

                if (n == 0)
                {
                    k1[0] = w[0];
                    k2[0] = -g/len*Math.Sin(theta[0]);
                    theta[1] = theta[0] + k1[0]*dt/2;
                    w[1]     = w[0]     + k2[0]*dt/2;
                }
                else if (n == 1)
                {
                    k1[1] = w[1];
                    k2[1] = -g/len*Math.Sin(theta[1]);
                    theta[2] = theta[0] + k1[1]*dt/2;
                    w[2]     = w[0]     + k2[1]*dt/2;
                }

                else if (n == 2)
                {
                    k1[2] = w[2];
                    k2[2] = -g/len*Math.Sin(theta[2]);
                    theta[3] = theta[0] + k1[2]*dt/2;
                    w[3]     = w[0]     + k2[2]*dt/2;
                }

                else
                {
                    k1[3] = w[3];
                    k2[3] = -g/len*Math.Sin(theta[3]);
                }

            }
            
            theta_answer = theta[0] + 1.0/6.0*(k1[0] + 2*k1[1] + 2*k1[2] + k1[3]) * dt;
            w_answer     = w[0]     + 1.0/6.0*(k2[0] + 2*k2[1] + 2*k2[2] + k2[3]) * dt;  

            theta[0] = theta_answer;
            w[0] = w_answer;
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