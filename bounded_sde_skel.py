import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.type_check import real
from scipy.linalg import toeplitz
from math import isclose


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
        solution_matrix = mm + column_vector + self.grav                    #############grav constant added
        solution_matrix = np.insert(solution_matrix, 0, self.v_0)                #The third parameter will be the initial value added, the second paramter is for position in array
        solution_matrix = solution_matrix.transpose()
       
        trap_matrix = (self.sigma * (1 + np.exp(self.b * self.deltat))  / 2   ) * upper_triangular
        matrix_multiplication = np.matmul(trap_matrix, random_vector )
        
        solution =  matrix_multiplication + column_vector - self.grav * self.deltat       #############grav constant added
        solution = np.insert(solution, 0, self.v_0)
        solution = solution.transpose()

        return solution

    def iteration(self):
        N = self.N
        x = np.zeros( (N + 1, 1) )
        solution = self.main_equation()
        dummy = 56 #random number that it should to actually hit
        for i in range(1, N+1): 
            x[i] = np.round(x[i - 1] + ( solution[i-1] * self.deltat), 4)

            tol = .00005

            if x[i] <= 0.1000 and x[i] >= -0.1000 :
                
                # dummy vairable is created in elif statements below but is used here in order to
                #skip an iteration that is created there....

                if dummy == x[i]:
                    continue
                else:
                    x[i] = x[i - 1] + ( solution[i] * self.deltat) 

            elif x[i] < -0.1000:
                h = 5000
                list = np.zeros((h+1, 1))
                time = np.linspace(solution[i-1], solution[i], h)
                z = int(h/2)
                for j in range(1, h):
                    list[j] = np.round(x[i-1] + (time[j])*(self.T/h)*j, 3)
                for g in range(1,z):
                    if list[z] < -0.1000:
                        #only thing changing is the list[z,g] indexing
                        list[z - g] = np.round(x[i-1] + (time[(z - g)])*(self.T/h)*(z - g), 3)
                        if isclose(list[z-g], -0.1000, rel_tol= tol):
                            x[i] = np.round(x[i-1] + (list[z-g] * (self.T/h)*(z - g)), 3)
                            x[i+1] = np.round(x[i] + (-list[z-g] * (self.deltat - (self.T/h)*(z - g))), 3)
                            self.v0 = -list[z-g]
                            solution = self.main_equation()
                            dummy = x[i+1]
                            
                    else:
                        list[z + g] = np.round(x[i-1] + (time[z + g])*(self.T/h)*(z + g), 3)
                        if isclose(list[z+g], -0.1000, rel_tol= tol):
                            x[i] = np.round(x[i-1] + (list[z+g] * (self.T/h)*(z + g)), 3)
                            x[i+1] = np.round(x[i] + (-list[z+g] * (self.deltat - (self.T/h)*(z + g))), 3)
                            self.v0 = -list[z+g]
                            solution = self.main_equation()
                            dummy = x[i+1]
                            
                                  

            elif x[i] > 0.1000:
                h = 5000
                list1 = np.zeros((h+1, 1))
                time1 = np.linspace(solution[i-1], solution[i], h)
                z = int(h/2)
                for j in range(1, h):
                    list1[j] = np.round(x[i-1] + (time1[j])*(self.T/h)*j, 3)
                for g in range(1,z):
                    if list1[z] < 0.1000:
                        list1[z + g] = np.round(x[i-1] + (time1[z+g])*(self.T/h)*(z+g), 3)
                        if isclose(list1[z+g], 0.1000, rel_tol= tol):
                            x[i] = np.round(x[i-1] + (list1[z+g] * (self.T/h)*(z+g)), 3)
                            x[i+1] = np.round(x[i] + (-list1[z+g] * (self.deltat - (self.T/h)*(z+g))), 3)
                            self.v0 = -list1[z+g]
                            solution = self.main_equation()
                            dummy = x[i+1]
                    else:
                        list1[z - g] = np.round(x[i-1] + (time1[z-g])*(self.T/h)*(z-g), 3)
                        if isclose(list1[z-g], 0.1000, rel_tol= tol):
                            x[i] = np.round(x[i-1] + (list1[z-g] * (self.T/h)*(z - g)), 3)
                            x[i+1] = np.round( x[i] + (-list1[z-g] * (self.deltat - (self.T/h)*(z - g))), 3)
                            self.v0 = -list1[z-g]
                            solution = self.main_equation()
                            dummy = x[i+1]
                            


                # for j in range(1, h):
                #     list_1[j] = np.round(x[i] + (time_1[j])*(self.T/h)*j, 4)
                #     print(list_1[j])
        
                # while j < h:
                #     h = h/2
                #     if list[h] < -.1000:
                #         start = h/2
                #         list[start + i] = np.round(x[i] + (time[j])*(self.T/h)*j, 4)
                #         j += 1
                #     else:
                #         if list[h] > -.1000:
                #             end = h/2
                #             list[end - i] = np.round(x[i] + (time[j])*(self.T/h)*j, 4)
                #             j += 1





                    # if isclose(list[j], -0.10, rel_tol=.0009):
                    #     print(list[j])
                    # else:
                    #     break

                  


            
            # elif x[i] < -0.10000:
            #     xx = np.zeros( (Nn + 1, 1) )
            #     z_i_l = np.linspace(solution[i-1], solution[i], Nn)
            #     dum = np.full(Nn+1, x[i-1])
            #     record_t = i*self.deltat
            #     for j in range(1,Nn):
            #         xx[j] = dum[j-1] + ( z_i_l[j] * (record_t - (self.T/Nn)*j)) 
            #         if isclose(xx[i], -0.10, rel_tol=.09):
            #             record_t = (self.T/Nn)*j                 #this record_t will have the delta t for that step
            #             x[i] = x[i-1] + ( xx[j] * record_t)
            #             self.v_0 = -z_i_l[j]
            #             solution = self.main_equation()
            #             x[i+1] = x[i] + (solution[i] *(self.deltat - record_t)) #it should now be followinng th egrid after this
            #             break             
            #         else:
            #             continue

            # elif x[i] > 0.10000:
            #     xxx = np.zeros( (Nn + 1, 1) )
            #     zz_i_l = np.linspace(solution[i-1], solution[i], Nn)
            #     ddum = np.full(Nn+1, x[i-1])
            #     rrecord_t = i*self.deltat
            #     for j in range(1,Nn):
            #         xxx[j] = ddum[j-1] + ( zz_i_l[j] * (rrecord_t - (self.T/Nn)*j)) 
            
            #         if isclose(xxx[i], 0.10, rel_tol=.009):
            #             rrecord_t = (self.T/Nn)*j                 #this record_t will have the delta t for that step
            #             x[i] = x[i-1] + (zz_i_l[j] * (rrecord_t + (self.T/Nn)*j))
            #             print(x[i])
            #             self.v_0 = -zz_i_l[j]
            #             solution = self.main_equation()
            #             x[i+1] = x[i] + (solution[i] * (self.deltat - rrecord_t)) #it should now be followinng th egrid after this
            #             dummy = x[i+1]
            #             break             
            #         else:
            #             continue

          
                #print('i-1:',solution[i-1],'i:',solution[i] )

        return x
'''
            elif x[i] > 0.10000:
                            
                x[i] = x[i - 1] + ( -solution[i-1] * self.deltat)
                self.v_0 = -solution[i-1]
                solution = self.main_equation()
                velocity = self.main_equation()
               
            # incorporate intersection...(look where it hit the boundary) then at that time make it reflect back and hit ht enext time grid
            #go from un-uniform grid to uniform grid for that time step
            #####This will be the test implimentation to check to see if this works####
            # elif x[i] < -0.10000:
            ###Working Code below###
            elif x[i] < -0.10000:
                
                x[i] = x[i - 1] + ( -solution[i-1] * self.deltat)
                self.v_0 = solution[i-1]
                solution = self.main_equation()
                #vel_vec = self.main_equation()
               
                            
            else:
                print('error')
'''

   
class Graph(Values):   
        def graph(self):

            meshpoints = np.linspace(self.t0, self.T, self.N + 1).transpose() 
            plt.rcParams['lines.linewidth'] = 2
            colors = ["coral","silver", "burlywood","lightgreen", "plum"]

            for i in range(5):

                plt.figure(1, figsize=(12,10))
                plt.plot(meshpoints, self.iteration() , label = "Path " + str(i+1), color=colors[i]) 
                plt.title('Trapezoidal Method - Bounded Position')
                plt.xlabel('t', size=14)
                plt.ylabel('Function',size=14)
            
        
            plt.legend()
            plt.show()


def main():
    N = (1000)
    numerical_method = Values(0, 1, N, 0)                             #instantiates the above class and methods
    numerical_method.main_equation()

    plotting_class = Graph(0, 1, N, 0)
    plotting_class.graph()

if __name__ == "__main__":
    main()                                         #since running program as main file