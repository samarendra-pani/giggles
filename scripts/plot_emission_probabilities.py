import matplotlib.pyplot as plt
import numpy as np

def generate_temperature_plot():
    # Parameters mirroring your C++ configuration
    k = 100
    temperatures = [10, 20, 30, 40]
    lower_bound = 1e-30
    
    # Material Design color palette
    colors = ['#2196F3', '#4CAF50', '#FF9800', '#E91E63']
    
    # Set overall aesthetic
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    for idx, T in enumerate(temperatures):
        exp_T = np.exp(T)
        
        # Translating the C++ mathematical normalization bounds
        alpha = (exp_T * lower_bound - 1.0) / (exp_T - 1.0)
        beta = (1.0 - lower_bound) / (exp_T - 1.0)
        
        probs = np.zeros(k + 1)
        probs[0] = lower_bound
        
        # Compute probabilities matching C++ loop logic
        for i in range(1, k):
            g = i / k
            probs[i] = alpha + (beta * np.exp(g * T))
            
        probs[k] = 1.0
        i_vals = np.arange(k + 1)/k
        
        # Plot on Linear Scale
        ax1.plot(i_vals, probs, color=colors[idx], linewidth=2.5, label=f'Temp = {T}')
        
        # Plot on Logarithmic Scale
        ax2.plot(i_vals, probs, color=colors[idx], linewidth=2.5, label=f'Temp = {T}')

    # Apply Material Design layout properties to both subplots
    for ax in [ax1, ax2]:
        ax.set_xlabel('g ---> [1 - min(1, f)]', fontsize=12, fontweight='500', color='#555555')
        ax.legend(frameon=False, fontsize=11)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    # Specifics for Linear axis
    ax1.set_ylabel('Emission Probability', fontsize=12, fontweight='500', color='#555555')
    ax1.set_title('Probability Cache (Linear Scale)', fontsize=14, fontweight='bold', pad=15, color='#333333')

    # Specifics for Log axis
    ax2.set_yscale('log')
    ax2.set_ylabel('Emission Probability (Log Scale)', fontsize=12, fontweight='500', color='#555555')
    ax2.set_title('Probability Cache (Log Scale)', fontsize=14, fontweight='bold', pad=15, color='#333333')
    # Set a custom bottom limit so the 1e-30 bounds don't visually squash the higher values
    ax2.set_ylim(bottom=1e-18, top=2.0)

    # Render and save
    plt.tight_layout()
    # plt.savefig('temperature_curves.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    generate_temperature_plot()