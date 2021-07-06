import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz


class Values(object):

    def __init__(self, t0, T, N, v_0):
        self.t0 = t0
        self.T = T
        self.N = N
        self.deltat = self.T / self.N
        self.sigma = .6
        self.b = 4
        self.grav = 0.2
        self.v_0 = v_0


    def main_equation(self):

        matrix = toeplitz(c=[1, *np.arange(2, self.N + 1)], r=np.arange(1,self.N+1))
        matrix_exp = np.exp(-self.b * matrix * self.deltat)

        upper_triangular = np.tril(matrix_exp) 
        column_vector = np.ones((self.N, 1))
        for i in range(self.N):
            column_vector[i] = np.exp( -self.b * (i+1) * self.deltat)
        column_vector *= self.v_0

        mean, sd = 0, np.sqrt(self.deltat)  
        random_values = np.random.normal(mean, sd, self.N) 
        random_vector = np.reshape(random_values, (self.N, 1))

        mm = self.sigma * np.matmul(upper_triangular, random_vector)        #Added grav here                                                            
        solution_matrix = mm + column_vector + self.grav                   
        solution_matrix = np.insert(solution_matrix, 0, self.v_0)           #The third parameter will be the initial value added, the second paramter is for position in array
        solution_matrix = solution_matrix.transpose()
       
        trap_matrix = (self.sigma * (1 + np.exp(self.b * self.deltat))  / 2   ) * upper_triangular
        matrix_multiplication = np.matmul(trap_matrix, random_vector )
        
        solution =  matrix_multiplication + column_vector - self.grav * self.deltat       #grav constant added
        solution = np.insert(solution, 0, self.v_0)
        solution = solution.transpose()

        return solution

    def iteration(self):
        N = self.N
        x = np.zeros( (N + 1, 1) )
        solution = self.main_equation()
  
        for i in range(1, N+1): 
            x[i] = np.round(x[i - 1] + ( solution[i-1] * self.deltat), 6)
            if x[i] <= 0.1000 and x[i] >= -0.1000 :
                x[i] = x[i - 1] + ( solution[i] * self.deltat) 

            elif x[i] < -0.1000:
                x[i] = np.round(x[i - 1] + ( solution[i-1] * self.deltat), 6)
                #Make the new dt
                dt =  np.abs(((-.1000 -x[i-1]) / solution[i-1]))
                #Redefine x[i] for when it is at the boundary
                x[i] = np.round(x[i - 1] + ( solution[i-1] * dt), 4)
                #Now define new Initial Values since reflecting back
                self.v0 = -solution[i-1]
                solution = self.main_equation()

            elif x[i] > 0.1000:
                x[i] = np.round(x[i - 1] + ( solution[i-1] * self.deltat), 6)
                dt =  np.abs(((.1000 - x[i-1]) / solution[i-1]))
                x[i] = np.round(x[i - 1] + ( solution[i-1] * dt), 4)
                #Step it at the new founded time t since that wil give you the point in wihcihc it crosses or grid size
                self.v0 = -solution[i-1]
                solution = self.main_equation()
                         
        return x
   
class Graph(Values):   
        def graph(self):

            meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose() 
            y = np.full(1001, -.1000)
            y1 = np.full(1001, .1000)
            plt.rcParams['lines.linewidth'] = 2
            colors = ["coral","silver", "burlywood","lightgreen", "plum"]

            for i in range(5):

                plt.figure(1, figsize=(12,10))
                plt.plot(meshpoints, self.iteration() , label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Bounded Position')
                plt.xlabel('Time', size=14)
                plt.ylabel('Position',size=14)
            
            plt.figure(1, figsize=(12,10))
            plt.plot(meshpoints, y, color = 'silver', linestyle = 'dashed', linewidth = 1.15)
            plt.plot(meshpoints, y1, color = 'silver', linestyle = 'dashed', linewidth = 1.15)        

            plt.legend()
            plt.show()

def main():
    N = (1000)
    numerical_method = Values(0, 1, N, 0)                             #instantiates the above class and methods
    numerical_method.main_equation()

    plotting_class = Graph(0, 1, N, 0)
    plotting_class.graph()

if __name__ == "__main__":
    main()                                                            #since running program as main file