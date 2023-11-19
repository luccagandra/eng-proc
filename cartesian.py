import matplotlib.pyplot as plt
import numpy as np

def plot_single_transfer(x, y, number, values,type, offset=0.1, single=True):
    """
    Plot a point in a 2D Cartesian space, display a straight line segment around it,
    plot a number slightly offset from the center, and add labels at the ends of the segments.

    Parameters:
    - x: x-coordinate of the point
    - y: y-coordinate of the point
    - number: The number to be plotted slightly offset from the center
    - offset_x: Offset in the x-direction for the number (default is 0.1)
    - offset_y: Offset in the y-direction for the number (default is 0.1)
    - values: Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
    """

    # Plotting the point
    plt.scatter(x, y, color='black', label='Point')

    # Displaying a straight line segment around the point
    segment_length = 1  # Length of the segment in both x and y directions

    if type == 'reta':
        plt.plot([x - segment_length/2, x + segment_length/2], [y, y], color='red')
        plt.plot([x, x], [y - segment_length/2, y + segment_length/2], color='blue')

        if values[2] == None:
            text_hot = ""
        else:
            text_hot = f'Q{values[0]} = {values[2]}'

        if values[3] == None:
            text_cold = ""
        else:
            text_cold = f'F{values[1]} = {values[3]}'

        # Adding labels 
        plt.text(x - segment_length/2, y + offset, text_hot, color='red', ha='right', va='center', fontsize=10)
        plt.text(x + offset, y + segment_length/2, text_cold, color='blue', ha='center', va='bottom', fontsize=10)

        plt.text(x + segment_length/2, y + offset, f'{values[4]}', color='red', ha='left', va='center', fontsize=10)
        plt.text(x + offset, y - segment_length/2, f'{values[5]}', color='blue', ha='center', va='top', fontsize=10)

    elif type == 'hot': # Usando água
        plt.plot([x - segment_length/4, x + segment_length/4], [y - segment_length/4, y + segment_length/4], color='blue')
        plt.plot([x - segment_length/2, x + segment_length/2], [y, y], color='red')
        
        plt.text(x - segment_length/4 - offset, y - segment_length/4 - offset, values[5], color='blue', ha='left', va='center', fontsize=10)
        plt.text(x + segment_length/4 + offset, y + segment_length/4 + offset, values[3], color='blue', ha='right', va='center', fontsize=10)
        plt.text(x + segment_length/2, y + offset, values[4], color='red', ha='left', va='center', fontsize=10)

    elif type == 'cold': # Usando vapor
        plt.plot([x - segment_length/4, x + segment_length/4], [y - segment_length/4, y + segment_length/4], color='red')
        plt.plot([x, x], [y - segment_length/2, y + segment_length/2], color='blue')

        plt.text(x - segment_length/4 - offset, y - segment_length/4 - offset, values[2], color='red', ha='left', va='center', fontsize=10)
        plt.text(x + segment_length/4 + offset, y + segment_length/4 + offset, values[4], color='red', ha='right', va='center', fontsize=10)
        plt.text(x + offset, y - segment_length/2, values[5], color='blue', ha='center', va='top', fontsize=10)

    # Plotting the number slightly offset from the center
    plt.text(x + offset/2, y + offset/2, str(number), color='black', ha='center', va='center', fontsize=12, fontweight='bold')

    if single == True:
        # Adding labels, title, legend, and grid for better visualization
        plt.title('Rede Final')
        plt.grid(False)
        plt.axis('off')
        plt.show()

def plot_multiple_transfers(plot_list, last_matrix):
    # Plotting the point
    plt.scatter(1, 1, color='red', label='Point')

    x,y = 0,0
    plotted_out = []
    values_to_clear = []

    for i in range(len(last_matrix)):
        if last_matrix[i][1] == last_matrix[i][2]:
            values_to_clear.append(last_matrix[i][1])

    for i in range(len(plot_list)):

        values = plot_list[i]

        if i == 0:
            plotted_out.append([values[4], 0,0, 'hot']) # TSQ, plotted at 0,0
            plotted_out.append([values[5], 0,0, 'cold']) # TSF, plotted at 0,0

            plot_single_transfer(0,0,(i+1),values,'reta', single=False)
            prev_value_x = 0
            prev_value_y = 0
            prev_values = values

        else:
            if values[3] == prev_values[5]: # Conexão por corrente fria
                x = prev_value_x
                y = prev_value_y - 1
                values_to_clear.append(values[3])
                values[3] = None

            if values[2] == prev_values[4]: # Conexão por corrente quente
                x = prev_value_x + 1
                y = prev_value_y
                values_to_clear.append(values[2])
                values[2] = None

            prev_value_x = x
            prev_value_y = y
            prev_values = values

            if i == len(plot_list)-1: # Última iteração
                plotted_out.append([values[4], x,y, 'hot']) 
                plotted_out.append([values[5], x,y, 'cold']) 

                plot_single_transfer(x,y,(i+1),values,'reta', single=False) # Terão elementos livres!

                #print(plotted_out)
                #print(values_to_clear)

                for j in range(len(plotted_out)-1, -1, -1):  
                    if plotted_out[j][0] in values_to_clear:
                        plotted_out.pop(j)
                
                for a in range(len(plotted_out)): # Completando com utilidades
                    
                    if a+1 == len(plotted_out):
                        single = True
                    else:
                        single = False  
                    
                    if plotted_out[a][3] == 'hot': # água de 30 a 50
                        
                        for k in range(len(last_matrix)):
                            if last_matrix[k][1] == plotted_out[a][0]:
                                meta = last_matrix[k][2]
                                
                        values = [0,0,plotted_out[a][0],30,meta,50] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
                        plot_single_transfer(plotted_out[a][1]+1, plotted_out[a][2], (a+i+2),values,'hot', single=single)

                    if plotted_out[a][3] == 'cold': # vapor de 250 a 250

                        for h in range(len(last_matrix)):
                            if last_matrix[h][1] == plotted_out[a][0]:
                                meta = last_matrix[h][2]

                        values = [0,0,250,plotted_out[a][0],250,meta] # Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
                        plot_single_transfer(plotted_out[a][1], plotted_out[a][2]-1, (a+i+2),values,'cold', single=single)

            else:
                plotted_out.append([values[4], x,y, 'hot']) 
                plotted_out.append([values[5], x,y, 'cold']) 
                plot_single_transfer(x,y,(i+1),values,'reta', single=False)
            


# Example usage:
#plot_single_transfer(1, 1, 1, offset=0.15)

# Qx, Fx, Qx_valor, Fx_valor, TSQ_valor, TSF_valor
plot_list = [[2, 2, 250.0, 100.0, 140.0, 131.4], [1, 2, 180.0, 131.4, 153.0, 170.0], [1, 1, 153.0, 60.0, 111.5, 143.0]]

last_matrix = [[10.0, 111.5, 90.0], [2.0, 140.0, 140.0], [5.0, 143.0, 150.0], [7.0, 170.0, 220.0]]

plot_multiple_transfers(plot_list, last_matrix)