import numpy as np
import matplotlib.pyplot as plt
from numpy.core.defchararray import endswith
from numpy.lib.twodim_base import diag
from scipy.linalg import toeplitz


class Values(object):

    def __init__(self, t0, T, N):
        self.t0 = t0
        self.T = T
        self.N = N
        self.deltat = self.T / self.N
        self.v_0 = 0
        self.sigma = .6
        self.b = 4
        self.grav = .2

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
        mm = self.sigma * np.matmul(self.uppertriangular(), self.delta_w_test())       #Added grav here                                                            
        solution_matrix = mm + self.exponential_matrix() + self.grav   #############grav constant added
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
        solution =  matrix_multiplication + self.exponential_matrix() + self.grav #############grav constant added
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
    #This will be using the trapezoidal method   
    def trapezoidal_iteration(self):
         x_trap = np.zeros( (self.N + 1, 1) )
         x_trap[0] = 0
         velocity_vector = self.matrix()
         for i in range(1, self.N +1):
             x_trap[i] = x_trap[i - 1] + ( (velocity_vector[i] + velocity_vector[i - 1]) / 2 ) * self.deltat
         return x_trap



class bounded_position(trapezoidal, Values):
    def left_point(self):
        x = np.zeros( (self.N + 1, 1) )
        velocity_vector = self.solution()

        # for i in range(1, self.N+1): 
        #     x[i] = x[i - 1] + ( velocity_vector[i-1] * self.deltat)
        #     if x[i] <= 0.10000 and x[i] >= -0.10000 :
        #         x[i] = x[i - 1] + ( velocity_vector[i-1] * self.deltat)              
        #     elif x[i] > 0.10000:
        #         x[i] = x[i - 1] + ( -velocity_vector[i-1] * self.deltat)
        #     elif x[i] < -0.10000:
        #         x[i] = x[i - 1] + ( -velocity_vector[i-1] * self.deltat)
        
        #     else:
        #         print('ERROR')

        bound = .10000
        for i in range(1, self.N+1): 
          velocity = (velocity_vector[i] + velocity_vector[i - 1]) / 2
          x[i] = x[i - 1] + ( velocity * self.deltat)
        
        while max(abs(x)) > bound:
            x[x<-bound] = -x[x<-bound] -2*bound
            x[x>bound] = -x[x>bound] + 2*bound


        return x


    def bounded_trap(self):
        x = np.zeros( (self.N + 1, 1) )
        velocity_vector = self.matrix()
        bound = .10000
        for i in range(1, self.N+1): 
          velocity = (velocity_vector[i] + velocity_vector[i - 1]) / 2
          x[i] = x[i - 1] + ( velocity * self.deltat)
        
        while max(abs(x)) > bound:
            x[x<-bound] = -x[x<-bound] -2*bound
            x[x>bound] = -x[x>bound] + 2*bound


        return x




    #Expected Value for these values from above is 0 since E[X_0] is 0 and so is E[V_0]

    #Expected Value and Variance 
    # def expected_leftpoint(self):
    #     x_0 = 0 #This is initial position
    #     meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose()
    #     list = np.zeros( (self.N + 1, 1) )

    #     for i in range(self.N + 1):
    #         list[i] = x_0 + ( ( 1 - np.exp(-self.b * (i*self.deltat)) ) / self.b ) * self.v_0
    #     return list


    # def variance_leftpoint(self):
    #     meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose()
    #     list = np.zeros( (self.N + 1, 1) )

    #     for i in range(self.N + 1):
    #         list[i] = ((self.sigma**2)*(i*self.deltat) / self.b**2) - ( self.sigma**2 / (2*self.b**3))*(3 - 4*np.exp(-self.b**(i*self.deltat)) + np.exp(-2*self.b*(i*self.deltat))) 
    #     return list



    
class Graph(bounded_position, position, trapezoidal, Values):                                                    
    def graph(self):
        #np.random.seed(0)
        meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose() 
        plt.rcParams['lines.linewidth'] = 1
        colors = ["coral","silver", "burlywood","lightgreen", "plum"]
        mean_matrix = self.exponential_matrix()
        mean_matrix = np.insert(mean_matrix, 0, 0) #Hard coding first expected value here
        for i in range(5):
            try:
                # plt.figure(1, figsize=(12,10))
                # plt.plot(meshpoints, self.solution(), label = "Path " + str(i+1), color=colors[i]) 
                # plt.title('Left-Point Method - Velocity')
                # plt.xlabel('t', size=14)
                # plt.ylabel('Function',size=14)
                # plt.legend()
                    
                # plt.figure(2, figsize=(12,10) )
                # plt.plot(meshpoints, self.matrix(), label = "Path " + str(i+1), color=colors[i])
                # plt.title('Trapezoidal Method - Velocity')
                # plt.xlabel('t', size=14)
                # plt.ylabel('Function',size=14)
                
                # # plt.figure(3, figsize=(12,10))
                # # plt.plot(meshpoints, self.lp_iteration(), label = "Path " + str(i+1), color=colors[i]) 
                # # plt.title('Left-Point Method - Position')
                # # plt.xlabel('t', size=14)
                # # plt.ylabel('Function',size=14)
                # # plt.legend()
                
                # plt.figure(3, figsize=(12,10))
                # plt.plot(meshpoints, self.trapezoidal_iteration(), label = "Path " + str(i+1), color=colors[i]) 
                # plt.title('Trapezoidal Method - Position')
                # plt.xlabel('t', size=14)
                # plt.ylabel('Function',size=14)
                # plt.legend()

                plt.figure(3, figsize=(12,10))
                plt.plot(meshpoints, self.left_point(),  label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Left Point - Bounded Position')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                plt.legend()

                plt.figure(4, figsize=(12,10))
                plt.plot(meshpoints, self.bounded_trap(), label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Trapezoidal Method - Bounded Position')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
                plt.legend()

            except ValueError:                                              
                print('~~~ARRAYS ARE NOT THE SAME SIZE~~~')
    
        ##                          ##
        # plt.figure(1, figsize=(12,10))
        # plt.plot(meshpoints, mean_matrix, color = 'grey', linestyle = 'dashed')
        # plt.plot(meshpoints, np.sqrt(self.variance()), color = 'grey', linestyle = 'dashed')
    
        # plt.figure(2, figsize=(12,10))
        # plt.plot(meshpoints, mean_matrix, color = 'grey', linestyle = 'dashed')
        # plt.plot(meshpoints, np.sqrt(self.variance_trapezoidal()), color = 'grey', linestyle = 'dashed')

        ##the above are still used##

        # plt.figure(3, figsize=(12,10))
        # plt.plot(meshpoints, self.expected_leftpoint(), color = 'grey', linestyle = 'dashed')
        # plt.plot(meshpoints, np.sqrt(self.variance_leftpoint()), color = 'grey', linestyle = 'dashed')
 
        plt.legend()
        