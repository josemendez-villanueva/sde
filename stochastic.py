import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

class Values(object):

    def __init__(self, t0, T, N):
        self.t0 = t0
        self.T = T
        self.N = N
        self.deltat = self.T / self.N
        self.v_0 = 1
        self.sigma = 1
        self.b = 1

    def Matrix(self):
        matrix = np.ones((self.N, self.N))
        i = 0                                                   # for loop dummy vairable to act as a condition/counter
        t = 0                                                   # while loop dummy variable to act as a condition/counter
        while i < self.N:                                       #This creates the Matrix
            for row in matrix:
                if i == t:
                    for j in range(self.N):
                            matrix[i][j] = (np.exp(-self.b* ( (i+1)+(-j) )* self.deltat))
                    i += 1        
                else:
                    break
            t += 1
        return matrix

    def uppertriangular(self):
        upper_triangular = np.tril(self.Matrix())               #This creates our Upper Triangular matrix, taking away the unnecessary values
        return upper_triangular
        

    def exponential_matrix(self):                               #This column vector is multiplied by our V_0 initial condition
        column_vector = np.ones((self.N, 1))
        for i in range(self.N):
            column_vector[i] = np.exp( -self.b * (i+1) * self.deltat)
        column_vector *= self.v_0
        return column_vector

    # def pdf(self):                                              #This is for the pdf function

    #     meshpoints = np.linspace(self.t0, self.T, self.N + 1)
    #     column_vector = np.zeros((self.N + 1, 1))
    #     norm.pdf(x=0, loc= 0, scale = self.deltat)                          #loc is mean, scale is variance
    #     j = 0
    #     for i in range(len(meshpoints)):   
    #         if i == j:
    #             column_vector[i] = norm.pdf(meshpoints[j])
    #             j += 1
    #         else:  
    #             continue
    #     #column_vector = column_vector[:-1]
    #     return column_vector
   

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


    # def variance(self):
    #     variance_vector = np.ones((self.N + 1, 1))
    #     for i in range(self.N + 1):
    #         variance_vector[i] = (self.sigma ** 2) * (((np.exp( (-2 * self.b * self.deltat) * ( 1 - np.exp(-2 * (i) * self.b * self.deltat ))))) / (1 - np.exp(-2 * self.b * self.deltat))) * self.deltat
    #     return variance_vector


class trapezoidal(Values):
    
    def matrix(self):
        trap_matrix = (self.sigma * (1 + np.exp(self.b * self.deltat))  / 2   ) * self.uppertriangular ()
        matrix_multiplication = np.matmul(trap_matrix, self.delta_w_test() )
        solution =  matrix_multiplication + self.exponential_matrix()
        solution = np.insert(solution, 0, self.v_0)
        solution = solution.transpose()

        return solution

    # def variance_trapezoidal(self):
    #     variance_vector = np.ones((self.N + 1, 1))
    #     for i in range(self.N + 1):
    #         variance_vector[i] = (self.sigma ** 2) * (((self.sigma**2) * (1 + np.exp(self.b * self.deltat)) / 4)) * (((np.exp( (-2 * self.b * self.deltat) * ( 1 - np.exp(-2 * (i) * self.b * self.deltat ))))) / (1 - np.exp(-2 * self.b * self.deltat))) * self.deltat

    #     return variance_vector
    



class Graph(trapezoidal, Values):                                                    #Created this class so when Numerical Method 2 is implemented it will inherit from it as well
    
    def graph(self):

        meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose() #This creats our equally spaced t values needed for our plots

        plt.rcParams['lines.linewidth'] = 3
        colors = ["coral","silver", "burlywood","lightgreen", "plum"]
        mean_matrix = self.exponential_matrix()
        mean_matrix = np.insert(mean_matrix, 0, np.exp(0))


        for i in range(5):
        

            try:
                plt.figure(1, figsize=(12,10))
                plt.plot(meshpoints, self.solution(),  label = "Path " + str(i+1), color=colors[i])
                plt.title('Left-Point Method')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                plt.legend()
                    
                plt.figure(2, figsize=(12,10) )
                plt.plot(meshpoints, self.matrix(), label = "Path " + str(i+1), color= colors[i])
                plt.title('Trapezoidal Method')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                
                          
            except ValueError:                                              #This tends to be the most common mistake for me...wrong sized arrays
                print('~~~ARRAYS ARE NOT THE SAME SIZE~~~')

        plt.figure(1, figsize=(12,10))
        plt.plot(meshpoints, mean_matrix, color = 'grey', linestyle = 'dashed')
        #plt.plot(meshpoints, self.variance(), color = 'grey', linestyle = 'dashed' ) 




        plt.figure(2, figsize=(12,10))
        plt.plot(meshpoints, mean_matrix, color = 'grey', linestyle = 'dashed')
        #plt.plot(meshpoints, self.variance_trapezoidal(), color = 'grey', linestyle = 'dashed' )

 
        plt.legend()
        plt.show()




def main():
    numerical_method = Values(0, 1, 500)                             #instantiates the above class and methods
    numerical_method.solution()
    trap = trapezoidal(0, 1, 500)
    trap.matrix()
    #trap.variance_trapezoidal()
    plotting_class = Graph(0, 1, 500)
    plotting_class.graph()

if __name__ == "__main__":                                         #since running program as main file
    main()
