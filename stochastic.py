import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import toeplitz





start = 0
end = 5
N = 500


class Values(object):

    def __init__(self, t0, T, N):
        self.t0 = t0
        self.T = T
        self.N = N
        self.deltat = self.T / self.N
        self.v_0 = 1
        self.sigma = .6
        self.b = 4

    def Matrix(self):
        
        matrix = toeplitz(c=[1, *np.arange(2, self.N + 1)], r=np.arange(1,self.N+1))
        matrix_exp = np.exp(-self.b * matrix * self.deltat)
        return matrix_exp

        
    def uppertriangular(self):
        upper_triangular = np.tril(self.Matrix())               #This creates our Upper Triangular matrix, taking away the unnecessary values
        return upper_triangular
        

    def exponential_matrix(self):                               #This column vector is multiplied by our V_0 initial condition
        column_vector = np.ones((self.N, 1))
        for i in range(self.N):
            column_vector[i] = np.exp( -self.b * (i+1) * self.deltat)
        column_vector *= self.v_0
        return column_vector


    def delta_w_test(self): 
        mean, sd = 0, np.sqrt(self.deltat) #following command uses sd and not variance so take sqrt
        random_values = np.random.normal(mean, sd, self.N)
        random_vector = np.reshape(random_values, (self.N, 1))
        return random_vector


    def solution(self):                                                                                                                       
        mm = self.sigma * np.matmul(self.uppertriangular(), self.delta_w_test())                                                                      
        solution_matrix = mm + self.exponential_matrix()
        solution_matrix = np.insert(solution_matrix, 0, self.v_0)   #The third parameter will be the initial value added, the second paramter is for position in array
        solution_matrix = solution_matrix.transpose()
        return solution_matrix


    def variance(self):
        variance_vector = np.ones((self.N + 1, 1))
        for i in range(self.N + 1):
            variance_vector[i] = (self.sigma ** 2) * (((np.exp (-2 * self.b * self.deltat) * ( 1 - np.exp(-2 * (i) * self.b * self.deltat )))) / (1 - np.exp(-2 * self.b * self.deltat))) * self.deltat
        return variance_vector


class trapezoidal(Values):
    
    def matrix(self):
        trap_matrix = (self.sigma * (1 + np.exp(self.b * self.deltat))  / 2   ) * self.uppertriangular ()
        matrix_multiplication = np.matmul(trap_matrix, self.delta_w_test() )
        solution =  matrix_multiplication + self.exponential_matrix()
        solution = np.insert(solution, 0, self.v_0)
        solution = solution.transpose()

        return solution

    def variance_trapezoidal(self):
        variance_vector = np.ones((self.N + 1, 1))
        for i in range(self.N + 1):
            variance_vector[i] = (self.sigma ** 2) * (( (1 + np.exp(self.b * self.deltat) ** 2) / 4)) * (((np.exp (-2 * self.b * self.deltat) * ( 1 - np.exp(-2 * (i) * self.b * self.deltat )))) / (1 - np.exp(-2 * self.b * self.deltat))) * self.deltat

        return variance_vector
    
        
class position(trapezoidal, Values):
    
        
    #This will be using the left point method
        
    def lp_iteration(self):
        velocity_vector = self.solution()
        
    
        x = np.zeros( (self.N + 1, 1) )
        
        
        x[0] = 0   #This is the initial position value
        

        for i in range(1, self.N + 1 ):
            x[i] =  x[i - 1] + (velocity_vector[i - 1]) * self.deltat 
        return x
        
    def trapezoidal_iteration(self):
        
        
         x_trap = np.zeros( (self.N + 1, 1) )
         x_trap[0] = 0
        
         velocity_vector = self.matrix()
         
         for i in range(1, self.N +1):
             x_trap[i] = x_trap[i - 1] + ( (velocity_vector[i] + velocity_vector[i - 1]) / 2 ) * self.deltat
    
         return x_trap
    
    
    
    

class Graph(position, trapezoidal, Values):                                                    
    
    def graph(self):

        meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose() 

        plt.rcParams['lines.linewidth'] = 3
        colors = ["coral","silver", "burlywood","lightgreen", "plum"]
        mean_matrix = self.exponential_matrix()
        mean_matrix = np.insert(mean_matrix, 0, np.exp(0))


        for i in range(5):
        

            try:
                plt.figure(1, figsize=(12,10))
                plt.plot(meshpoints, self.solution(), label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Left-Point Method - Velocity')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                plt.legend()
                    
                plt.figure(2, figsize=(12,10) )
                plt.plot(meshpoints, self.matrix(), label = "Path " + str(i+1), color=colors[i])
                plt.title('Trapezoidal Method - Velocity')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                
                plt.figure(3, figsize=(12,10))
                plt.plot(meshpoints, self.lp_iteration(), label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Left-Point Method - Position')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                plt.legend()
                
                plt.figure(4, figsize=(12,10))
                plt.plot(meshpoints, self.trapezoidal_iteration(), label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Trapezoidal Method - Position')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                plt.legend()
                
        
            except ValueError:                                              
                print('~~~ARRAYS ARE NOT THE SAME SIZE~~~')

        plt.figure(1, figsize=(12,10))
        plt.plot(meshpoints, mean_matrix, color = 'grey', linestyle = 'dashed')
        plt.plot(meshpoints, np.sqrt(self.variance()), color = 'grey', linestyle = 'dashed')
    
        plt.figure(2, figsize=(12,10))
        plt.plot(meshpoints, mean_matrix, color = 'grey', linestyle = 'dashed')
        plt.plot(meshpoints, np.sqrt(self.variance_trapezoidal()), color = 'grey', linestyle = 'dashed')
 
        plt.legend()
        plt.show()



def main():
    numerical_method = Values(start, end, N)                            
    numerical_method.solution()
    
    trap = trapezoidal(start, end, N)
    trap.matrix()
    
    X = position(start,end, N)
    X.lp_iteration()
    
    plotting_class = Graph(start, end, N)
    plotting_class.graph()

if __name__ == "__main__":                                        
    main()