#!/usr/bin/env python3
"""
Charge Analysis and Visualization Tools
========================================

Tools for analyzing and visualizing partial charges from Espaloma or other ML models.

Features:
- 3D charge distribution visualization
- Electrostatic potential surfaces
- Charge flow animations during vibrations
- Reactivity site prediction
- pKa estimation
- Hydrogen bonding analysis
- Charge transfer analysis
- Dipole moment decomposition

Usage:
    from charge_analysis import ChargeAnalyzer
    
    analyzer = ChargeAnalyzer(atoms, charges)
    analyzer.visualize_charges()
    analyzer.predict_reactive_sites()
    analyzer.analyze_all()

Author: Your Name
Date: 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from scipy.spatial.distance import cdist
from ase.neighborlist import natural_cutoffs, NeighborList
from typing import Optional, Tuple, List
import warnings

warnings.filterwarnings('ignore')

try:
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    pass


class ChargeAnalyzer:
    """
    Comprehensive charge analysis and visualization toolkit.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object containing molecular structure
    charges : np.ndarray
        Partial charges for each atom (in elementary charge units)
    """
    
    def __init__(self, atoms, charges: np.ndarray):
        """
        Initialize analyzer with molecular structure and charges.
        """
        self.atoms = atoms
        self.charges = charges
        self.symbols = atoms.get_chemical_symbols()
        self.positions = atoms.get_positions()
        self.natoms = len(atoms)
        
        # Build neighbor list for bond analysis
        self._build_neighbor_list()
        
        # Validate inputs
        if len(charges) != self.natoms:
            raise ValueError(f"Number of charges ({len(charges)}) doesn't match "
                           f"number of atoms ({self.natoms})")
    
    def _build_neighbor_list(self):
        """Build neighbor list for bonding analysis."""
        cutoffs = natural_cutoffs(self.atoms)
        self.nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        self.nl.update(self.atoms)
    
    # =========================================================================
    # VISUALIZATION FUNCTIONS
    # =========================================================================
    
    def visualize_charges(self, filename: str = 'charge_distribution.png',
                         show: bool = True, dpi: int = 300):
        """
        Create 3D visualization of charge distribution.
        
        Parameters
        ----------
        filename : str
            Output filename for saving the plot
        show : bool
            Whether to display the plot
        dpi : int
            Resolution for saved image
        """
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Color map: red = positive, blue = negative
        colors = ['red' if q > 0 else 'blue' for q in self.charges]
        sizes = [abs(q) * 1000 for q in self.charges]  # Size proportional to |charge|
        
        # Plot atoms with color-coded charges
        scatter = ax.scatter(
            self.positions[:, 0], 
            self.positions[:, 1], 
            self.positions[:, 2],
            c=self.charges, 
            cmap='RdBu_r', 
            s=sizes, 
            alpha=0.7,
            edgecolors='black', 
            linewidths=2,
            vmin=-max(abs(self.charges)), 
            vmax=max(abs(self.charges))
        )
        
        # Add labels for each atom
        for i, (pos, sym, q) in enumerate(zip(self.positions, self.symbols, self.charges)):
            ax.text(pos[0], pos[1], pos[2], 
                   f'{sym}{i}\n{q:+.2f}e', 
                   fontsize=9, ha='center', va='center',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))
        
        # Draw bonds
        self._draw_bonds(ax)
        
        # Color bar
        cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
        cbar.set_label('Partial Charge (e)', fontsize=12, fontweight='bold')
        
        # Labels and title
        ax.set_xlabel('X (x)', fontsize=12, fontweight='bold')
        ax.set_ylabel('Y (y)', fontsize=12, fontweight='bold')
        ax.set_zlabel('Z (z)', fontsize=12, fontweight='bold')
        ax.set_title('Atomic Partial Charges', fontsize=16, fontweight='bold', pad=20)
        
        # Set aspect ratio
        max_range = np.array([
            self.positions[:, 0].max() - self.positions[:, 0].min(),
            self.positions[:, 1].max() - self.positions[:, 1].min(),
            self.positions[:, 2].max() - self.positions[:, 2].min()
        ]).max() / 2.0
        
        mid_x = (self.positions[:, 0].max() + self.positions[:, 0].min()) * 0.5
        mid_y = (self.positions[:, 1].max() + self.positions[:, 1].min()) * 0.5
        mid_z = (self.positions[:, 2].max() + self.positions[:, 2].min()) * 0.5
        
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        
        plt.tight_layout()
        plt.savefig(filename, dpi=dpi, bbox_inches='tight')
        print(f"\u2713 Charge distribution plot saved to {filename}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def _draw_bonds(self, ax):
        """Draw bonds between atoms."""
        for i in range(self.natoms):
            indices, offsets = self.nl.get_neighbors(i)
            for j in indices:
                if i < j:  # Avoid drawing bonds twice
                    pos1 = self.positions[i]
                    pos2 = self.positions[j]
                    ax.plot([pos1[0], pos2[0]], 
                           [pos1[1], pos2[1]], 
                           [pos1[2], pos2[2]], 
                           'k-', linewidth=2, alpha=0.3)
    
    def plot_esp_surface(self, filename: str = 'esp_surface.html',
                        grid_points: int = 30):
        """
        Plot molecular electrostatic potential (MEP) surface.
        
        Parameters
        ----------
        filename : str
            Output filename (HTML for interactive plot)
        grid_points : int
            Number of grid points per dimension (larger = finer, slower)
        """
        try:
            import plotly.graph_objects as go
        except ImportError:
            print("\u26a0 Plotly not installed. Install with: pip install plotly")
            print("  Skipping ESP surface plot.")
            return
        
        print("Calculating electrostatic potential surface...")
        
        # Create 3D grid around molecule
        padding = 5  # Angstroms
        min_coords = self.positions.min(axis=0) - padding
        max_coords = self.positions.max(axis=0) + padding
        
        x = np.linspace(min_coords[0], max_coords[0], grid_points)
        y = np.linspace(min_coords[1], max_coords[1], grid_points)
        z = np.linspace(min_coords[2], max_coords[2], grid_points)
        X, Y, Z = np.meshgrid(x, y, z)
        
        # Calculate ESP at each grid point: V(r) = \u03a3\u1d62 q\u1d62/|r - r\u1d62|
        grid_points_array = np.c_[X.ravel(), Y.ravel(), Z.ravel()]
        distances = cdist(grid_points_array, self.positions)
        
        # Avoid singularities at atomic positions
        distances = np.maximum(distances, 0.1)
        
        # Calculate potential
        esp = np.sum(self.charges / distances, axis=1)
        ESP = esp.reshape(X.shape)
        
        # Create isosurface plot
        fig = go.Figure(data=go.Isosurface(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            value=ESP.flatten(),
            isomin=-0.15,
            isomax=0.15,
            surface_count=15,
            colorscale='RdBu_r',
            caps=dict(x_show=False, y_show=False, z_show=False),
            colorbar=dict(
                title=dict(text='ESP (e/Å)', side='right')  # ← FIXED
            )
        ))
        
        # Add atomic positions
        fig.add_trace(go.Scatter3d(
            x=self.positions[:, 0],
            y=self.positions[:, 1],
            z=self.positions[:, 2],
            mode='markers+text',
            marker=dict(size=8, color='black'),
            text=self.symbols,
            textposition='top center',
            name='Atoms'
        ))
        
        fig.update_layout(
            title='Electrostatic Potential Surface',
            scene=dict(
                xaxis_title='X (x)',
                yaxis_title='Y (y)',
                zaxis_title='Z (z)',
                aspectmode='data'
            ),
            width=900,
            height=700
        )
        
        fig.write_html(filename)
        print(f"\u2713 ESP surface saved to {filename} (open in web browser)")
    
    def animate_charge_flow(self, normal_mode: np.ndarray,
                           charge_calculator,
                           filename: str = 'charge_flow.gif',
                           frames: int = 50,
                           amplitude: float = 0.5):
        """
        Animate charge redistribution during molecular vibration.
        
        Parameters
        ----------
        normal_mode : np.ndarray
            Normal mode vector (3N array)
        charge_calculator : callable
            Function that calculates charges given an atoms object
        filename : str
            Output filename for animation
        frames : int
            Number of animation frames
        amplitude : float
            Maximum displacement amplitude in Angstroms
        """
        print(f"Creating charge flow animation ({frames} frames)...")
        
        equilibrium_pos = self.positions.copy()
        normal_mode_reshaped = normal_mode.reshape(-1, 3)
        
        # Normalize normal mode
        normal_mode_reshaped = normal_mode_reshaped / np.linalg.norm(normal_mode_reshaped)
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Animation frames: oscillate from -amplitude to +amplitude
        displacements = amplitude * np.sin(np.linspace(0, 2*np.pi, frames))
        
        def update(frame):
            ax.clear()
            
            # Displace along normal mode
            displaced_pos = equilibrium_pos + displacements[frame] * normal_mode_reshaped
            atoms_temp = self.atoms.copy()
            atoms_temp.set_positions(displaced_pos)
            
            # Calculate charges at this geometry
            try:
                _, charges_temp = charge_calculator.calculate_dipole(atoms_temp)
                if charges_temp is None:
                    charges_temp = self.charges  # Fallback
            except:
                charges_temp = self.charges  # Fallback if calculation fails
            
            # Plot
            sizes = [abs(q) * 1000 for q in charges_temp]
            scatter = ax.scatter(
                displaced_pos[:, 0], 
                displaced_pos[:, 1], 
                displaced_pos[:, 2],
                c=charges_temp, 
                cmap='RdBu_r', 
                s=sizes, 
                alpha=0.7,
                vmin=-max(abs(self.charges)), 
                vmax=max(abs(self.charges)),
                edgecolors='black',
                linewidths=1.5
            )
            
            ax.set_title(f'Charge Flow - Displacement: {displacements[frame]:.3f} �',
                        fontsize=14, fontweight='bold')
            ax.set_xlabel('X (x)', fontsize=11)
            ax.set_ylabel('Y (y)', fontsize=11)
            ax.set_zlabel('Z (z)', fontsize=11)
            
            # Keep axes fixed
            max_range = 5.0
            mid = equilibrium_pos.mean(axis=0)
            ax.set_xlim(mid[0] - max_range, mid[0] + max_range)
            ax.set_ylim(mid[1] - max_range, mid[1] + max_range)
            ax.set_zlim(mid[2] - max_range, mid[2] + max_range)
        
        anim = animation.FuncAnimation(fig, update, frames=frames, interval=50)
        
        try:
            anim.save(filename, writer='pillow', fps=20)
            print(f"\u2713 Charge flow animation saved to {filename}")
        except Exception as e:
            print(f"\u26a0 Could not save animation: {e}")
            print("  Make sure Pillow is installed: pip install Pillow")
        
        plt.close()
    
    # =========================================================================
    # ANALYSIS FUNCTIONS
    # =========================================================================
    
    def predict_reactive_sites(self, print_output: bool = True) -> dict:
        """
        Identify nucleophilic and electrophilic sites.
        
        Parameters
        ----------
        print_output : bool
            Whether to print results to console
        
        Returns
        -------
        dict
            Dictionary with keys 'nucleophiles', 'electrophiles', 'polarized_bonds'
        """
        nucleophiles = []
        electrophiles = []
        polarized_bonds = []
        
        # Identify nucleophilic sites (negative charges)
        for i, (sym, q) in enumerate(zip(self.symbols, self.charges)):
            if q < -0.3:
                nucleophiles.append({'index': i, 'symbol': sym, 'charge': q})
        
        # Identify electrophilic sites (positive charges, excluding H)
        for i, (sym, q) in enumerate(zip(self.symbols, self.charges)):
            if q > +0.3 and sym not in ['H']:
                electrophiles.append({'index': i, 'symbol': sym, 'charge': q})
        
        # Identify polarized bonds
        for i in range(self.natoms):
            indices, offsets = self.nl.get_neighbors(i)
            for j in indices:
                if i < j:  # Avoid double counting
                    q_diff = abs(self.charges[i] - self.charges[j])
                    if q_diff > 0.3:
                        polarized_bonds.append({
                            'atom1': i,
                            'atom2': j,
                            'symbol1': self.symbols[i],
                            'symbol2': self.symbols[j],
                            'charge1': self.charges[i],
                            'charge2': self.charges[j],
                            'delta_q': q_diff
                        })
        
        if print_output:
            print("\n" + "="*70)
            print("REACTIVE SITE ANALYSIS")
            print("="*70)
            
            print("\nNUCLEOPHILIC SITES (electron-rich, attacks \u03b4+):")
            if nucleophiles:
                for site in sorted(nucleophiles, key=lambda x: x['charge']):
                    print(f"  Atom {site['index']:2d} ({site['symbol']:2s}): "
                          f"{site['charge']:+.3f} e")
            else:
                print("  None found")
            
            print("\nELECTROPHILIC SITES (electron-poor, attacked by \u03b4-):")
            if electrophiles:
                for site in sorted(electrophiles, key=lambda x: x['charge'], reverse=True):
                    print(f"  Atom {site['index']:2d} ({site['symbol']:2s}): "
                          f"{site['charge']:+.3f} e")
            else:
                print("  None found")
            
            print("\nPOLARIZED BONDS (potential for heterolytic cleavage):")
            if polarized_bonds:
                for bond in sorted(polarized_bonds, key=lambda x: x['delta_q'], reverse=True):
                    print(f"  {bond['symbol1']}{bond['atom1']}-{bond['symbol2']}{bond['atom2']}: "
                          f"\u0394q = {bond['delta_q']:.3f} e "
                          f"({bond['symbol1']}: {bond['charge1']:+.2f}, "
                          f"{bond['symbol2']}: {bond['charge2']:+.2f})")
            else:
                print("  None found")
            
            print("="*70 + "\n")
        
        return {
            'nucleophiles': nucleophiles,
            'electrophiles': electrophiles,
            'polarized_bonds': polarized_bonds
        }
    
    def estimate_pKa(self, print_output: bool = True) -> List[dict]:
        """
        Estimate pKa values for acidic protons.
        
        Parameters
        ----------
        print_output : bool
            Whether to print results to console
        
        Returns
        -------
        list
            List of dictionaries with pKa estimates
        """
        pKa_estimates = []
        
        for i, (sym, q) in enumerate(zip(self.symbols, self.charges)):
            if sym == 'H':
                indices, _ = self.nl.get_neighbors(i)
                for j in indices:
                    neighbor_sym = self.symbols[j]
                    neighbor_q = self.charges[j]
                    
                    if neighbor_sym in ['O', 'N', 'S']:
                        # Empirical correlation (rough approximation!)
                        # More positive H = more acidic
                        # More negative heteroatom = stabilizes conjugate base
                        estimated_pKa = 15 - 20*q - 10*neighbor_q
                        
                        if estimated_pKa < 5:
                            acidity = "very acidic"
                        elif estimated_pKa < 10:
                            acidity = "moderately acidic"
                        elif estimated_pKa < 15:
                            acidity = "weakly acidic"
                        else:
                            acidity = "not acidic"
                        
                        pKa_estimates.append({
                            'H_index': i,
                            'X_index': j,
                            'X_symbol': neighbor_sym,
                            'H_charge': q,
                            'X_charge': neighbor_q,
                            'pKa': estimated_pKa,
                            'acidity': acidity
                        })
        
        if print_output:
            print("\n" + "="*70)
            print("ACIDIC PROTON ANALYSIS (rough pKa estimation)")
            print("="*70)
            print("\nNote: These are rough estimates based on partial charges.")
            print("      For accurate pKa, use specialized software or experimental data.")
            
            if pKa_estimates:
                for est in pKa_estimates:
                    print(f"\n{est['X_symbol']}{est['X_index']}-H{est['H_index']}:")
                    print(f"  H charge: {est['H_charge']:+.3f} e")
                    print(f"  {est['X_symbol']} charge: {est['X_charge']:+.3f} e")
                    print(f"  Estimated pKa: ~{est['pKa']:.1f} ({est['acidity']})")
            else:
                print("\nNo acidic protons (O-H, N-H, S-H) found")
            
            print("="*70 + "\n")
        
        return pKa_estimates
    
    def analyze_hydrogen_bonds(self, print_output: bool = True) -> dict:
        """
        Identify potential hydrogen bond donors and acceptors.
        
        Parameters
        ----------
        print_output : bool
            Whether to print results to console
        
        Returns
        -------
        dict
            Dictionary with keys 'donors' and 'acceptors'
        """
        donors = []
        acceptors = []
        
        # H-bond donors (H attached to O, N, F with positive charge)
        for i, (sym, q) in enumerate(zip(self.symbols, self.charges)):
            if sym == 'H' and q > 0.25:
                indices, _ = self.nl.get_neighbors(i)
                for j in indices:
                    if self.symbols[j] in ['O', 'N', 'F']:
                        donor_strength = q * abs(self.charges[j])
                        donors.append({
                            'H_index': i,
                            'X_index': j,
                            'X_symbol': self.symbols[j],
                            'H_charge': q,
                            'X_charge': self.charges[j],
                            'strength': donor_strength
                        })
        
        # H-bond acceptors (O, N, F with negative charge)
        for i, (sym, q) in enumerate(zip(self.symbols, self.charges)):
            if sym in ['O', 'N', 'F'] and q < -0.3:
                acceptor_strength = abs(q)
                acceptors.append({
                    'index': i,
                    'symbol': sym,
                    'charge': q,
                    'strength': acceptor_strength
                })
        
        if print_output:
            print("\n" + "="*70)
            print("HYDROGEN BONDING POTENTIAL")
            print("="*70)
            
            print("\nH-BOND DONORS:")
            if donors:
                for donor in sorted(donors, key=lambda x: x['strength'], reverse=True):
                    strength_desc = 'strong' if donor['strength'] > 0.2 else 'moderate'
                    print(f"  {donor['X_symbol']}{donor['X_index']}-H{donor['H_index']}: "
                          f"strength \u2248 {donor['strength']:.3f} ({strength_desc})")
            else:
                print("  None found")
            
            print("\nH-BOND ACCEPTORS:")
            if acceptors:
                for acceptor in sorted(acceptors, key=lambda x: x['strength'], reverse=True):
                    strength_desc = 'strong' if acceptor['strength'] > 0.5 else 'moderate'
                    print(f"  {acceptor['symbol']}{acceptor['index']}: "
                          f"strength \u2248 {acceptor['strength']:.3f} ({strength_desc})")
            else:
                print("  None found")
            
            print("="*70 + "\n")
        
        return {'donors': donors, 'acceptors': acceptors}
    
    def analyze_charge_transfer(self, groups: Optional[dict] = None,
                                print_output: bool = True) -> dict:
        """
        Quantify charge transfer between functional groups.
        
        Parameters
        ----------
        groups : dict, optional
            Dictionary mapping group names to lists of atom indices
            Example: {'methyl': [0, 1, 2, 3], 'carboxyl': [4, 5, 6, 7]}
        print_output : bool
            Whether to print results to console
        
        Returns
        -------
        dict
            Dictionary with group charges
        """
        total_charge = np.sum(self.charges)
        
        results = {'total': total_charge}
        
        if print_output:
            print("\n" + "="*70)
            print("CHARGE TRANSFER ANALYSIS")
            print("="*70)
            
            print(f"\nTotal molecular charge: {total_charge:+.6f} e")
            if abs(total_charge) < 0.01:
                print("(Close to zero - good, neutral molecule)")
            else:
                print("(Should match formal charge; deviation indicates numerical error)")
        
        if groups is not None:
            if print_output:
                print("\nGroup charges:")
            
            for group_name, indices in groups.items():
                group_charge = np.sum(self.charges[indices])
                results[group_name] = group_charge
                
                if print_output:
                    print(f"  {group_name}: {group_charge:+.3f} e")
            
            # Analyze transfer
            if len(groups) == 2 and print_output:
                group_names = list(groups.keys())
                transfer = -results[group_names[0]]
                print(f"\nCharge transfer from {group_names[0]} to {group_names[1]}: "
                      f"{transfer:.3f} e")
                if transfer > 0:
                    print(f"({group_names[0]} donates electrons to {group_names[1]})")
                else:
                    print(f"({group_names[1]} donates electrons to {group_names[0]})")
        
        if print_output:
            print("="*70 + "\n")
        
        return results
    
    def decompose_dipole(self, print_output: bool = True) -> dict:
        """
        Decompose dipole moment into atomic contributions.
        
        Parameters
        ----------
        print_output : bool
            Whether to print results to console
        
        Returns
        -------
        dict
            Dictionary with dipole information
        """
        # Calculate dipole moment
        ANGSTROM_TO_BOHR = 1.8897259886
        dipole_au = np.dot(self.charges, self.positions) / ANGSTROM_TO_BOHR  # e�Bohr
        dipole_magnitude = np.linalg.norm(dipole_au)
        dipole_debye = dipole_magnitude * 2.5418  # Convert to Debye
        
        # Calculate contributions
        contributions = np.array([q * pos / ANGSTROM_TO_BOHR 
                                 for q, pos in zip(self.charges, self.positions)])
        contrib_mags = np.linalg.norm(contributions, axis=1)
        
        if print_output:
            print("\n" + "="*70)
            print("DIPOLE MOMENT DECOMPOSITION")
            print("="*70)
            
            print(f"\nTotal dipole moment:")
            print(f"  Vector: [{dipole_au[0]:+.3f}, {dipole_au[1]:+.3f}, "
                  f"{dipole_au[2]:+.3f}] e�Bohr")
            print(f"  Magnitude: {dipole_magnitude:.3f} e�Bohr = {dipole_debye:.3f} Debye")
            
            if dipole_magnitude > 0.01:
                theta = np.arctan2(dipole_au[1], dipole_au[0]) * 180 / np.pi
                phi = np.arccos(dipole_au[2] / dipole_magnitude) * 180 / np.pi
                print(f"  Direction: \u03b8={theta:.1f}�, \u03c6={phi:.1f}�")
            
            print("\nAtomic contributions:")
            print(f"{'Atom':4s} {'Sym':3s} {'Charge':>7s} {'Position (�)':>30s} "
                  f"{'Contribution (e�Bohr)':>30s}")
            print("-" * 80)
            
            for i, (sym, q, pos, contrib) in enumerate(
                zip(self.symbols, self.charges, self.positions, contributions)):
                print(f"{i:3d}  {sym:2s}  {q:+.3f}   "
                      f"({pos[0]:6.2f}, {pos[1]:6.2f}, {pos[2]:6.2f})   "
                      f"({contrib[0]:+6.3f}, {contrib[1]:+6.3f}, {contrib[2]:+6.3f})")
            
            print("\nTop 5 contributors by magnitude:")
            top_contributors = np.argsort(contrib_mags)[::-1][:5]
            for rank, i in enumerate(top_contributors, 1):
                percent = (contrib_mags[i] / dipole_magnitude * 100) if dipole_magnitude > 0 else 0
                print(f"  {rank}. Atom {i} ({self.symbols[i]}): "
                      f"{contrib_mags[i]:.3f} e�Bohr ({percent:.1f}% of total)")
            
            print("="*70 + "\n")
        
        return {
            'dipole_vector': dipole_au,
            'dipole_magnitude': dipole_magnitude,
            'dipole_debye': dipole_debye,
            'contributions': contributions,
            'contribution_magnitudes': contrib_mags
        }
    
    # =========================================================================
    # COMPREHENSIVE ANALYSIS
    # =========================================================================
    
    def analyze_all(self, groups: Optional[dict] = None,
                   save_plots: bool = True,
                   output_dir: str = '.'):
        """
        Run all analysis and visualization functions.
        
        Parameters
        ----------
        groups : dict, optional
            Group definitions for charge transfer analysis
        save_plots : bool
            Whether to save visualization plots
        output_dir : str
            Directory for output files
        """
        import os
        
        if save_plots and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        print("\n" + "="*70)
        print(" "*20 + "COMPREHENSIVE CHARGE ANALYSIS")
        print("="*70 + "\n")
        
        # Analysis functions
        self.decompose_dipole(print_output=True)
        self.predict_reactive_sites(print_output=True)
        self.estimate_pKa(print_output=True)
        self.analyze_hydrogen_bonds(print_output=True)
        self.analyze_charge_transfer(groups=groups, print_output=True)
        
        # Visualizations
        if save_plots:
            print("Generating visualizations...")
            
            charge_file = os.path.join(output_dir, 'charge_distribution.png')
            self.visualize_charges(filename=charge_file, show=False)
            
            esp_file = os.path.join(output_dir, 'esp_surface.html')
            self.plot_esp_surface(filename=esp_file, grid_points=30)
            
            print("\n" + "="*70)
            print("ANALYSIS COMPLETE")
            print(f"Visualizations saved to: {output_dir}")
            print("="*70 + "\n")
    
    def analyze_all_to_file(self, filename: str = 'charge_analysis_report.txt',
                           groups: Optional[dict] = None,
                           save_plots: bool = True,
                           output_dir: str = '.'):
        """
        Run all analysis and save results to text file.
        
        Parameters
        ----------
        filename : str
            Output text file for analysis results
        groups : dict, optional
            Group definitions for charge transfer analysis
        save_plots : bool
            Whether to save visualization plots
        output_dir : str
            Directory for output files (plots and report)
        """
        import os
        import sys
        from io import StringIO
        
        if save_plots and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Create full path for report
        report_path = os.path.join(output_dir, filename)
        
        # Open file for writing
        with open(report_path, 'w', encoding='utf-8') as f:
            # Write header
            f.write("="*70 + "\n")
            f.write(" "*15 + "COMPREHENSIVE CHARGE ANALYSIS REPORT\n")
            f.write("="*70 + "\n")
            f.write(f"\nMolecule: {len(self.atoms)} atoms\n")
            f.write(f"Symbols: {' '.join(self.symbols)}\n")
            f.write(f"Total charge: {np.sum(self.charges):+.6f} e\n")
            f.write("\n" + "="*70 + "\n\n")
            
            # Redirect stdout to capture print statements
            old_stdout = sys.stdout
            sys.stdout = f
            
            try:
                # Run all analyses (output goes to file)
                print("SECTION 1: DIPOLE MOMENT ANALYSIS")
                print("="*70)
                self.decompose_dipole(print_output=True)
                
                print("\n" + "="*70)
                print("SECTION 2: REACTIVE SITES")
                print("="*70)
                self.predict_reactive_sites(print_output=True)
                
                print("\n" + "="*70)
                print("SECTION 3: ACIDITY ANALYSIS")
                print("="*70)
                self.estimate_pKa(print_output=True)
                
                print("\n" + "="*70)
                print("SECTION 4: HYDROGEN BONDING")
                print("="*70)
                self.analyze_hydrogen_bonds(print_output=True)
                
                print("\n" + "="*70)
                print("SECTION 5: CHARGE TRANSFER")
                print("="*70)
                self.analyze_charge_transfer(groups=groups, print_output=True)
                
                print("\n" + "="*70)
                print("END OF REPORT")
                print("="*70)
                
            finally:
                # Restore stdout
                sys.stdout = old_stdout
        
        # Print progress to console (not to file)
        print(f"  ✓ Analysis report saved to: {report_path}")
        
        # Generate visualizations
        if save_plots:
            print(f"  → Generating visualizations...")
            
            charge_file = os.path.join(output_dir, 'charge_distribution.png')
            self.visualize_charges(filename=charge_file, show=False, dpi=300)
            
            esp_file = os.path.join(output_dir, 'esp_surface.html')
            try:
                self.plot_esp_surface(filename=esp_file, grid_points=30)
            except ImportError:
                print(f"  ⚠ Skipping ESP surface (plotly not installed)")
            
            print(f"  ✓ Visualizations saved to: {output_dir}/")


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    """
    Example usage of ChargeAnalyzer
    """
    from ase.io import read
    
    print("\n" + "="*70)
    print("CHARGE ANALYSIS TOOLKIT - EXAMPLE USAGE")
    print("="*70 + "\n")
    
    # Example 1: Load molecule and charges
    print("Example: To use this toolkit with your gm_main.py code:\n")
    
    example_code = """
    from charge_analysis import ChargeAnalyzer
    from ase.io import read
    
    # After calculating charges with Espaloma:
    mol = read('molecule.xyz')
    dipole, charges = espaloma_calc.calculate_dipole(mol)
    
    # Create analyzer
    analyzer = ChargeAnalyzer(mol, charges)
    
    # Run individual analyses:
    analyzer.visualize_charges()
    analyzer.predict_reactive_sites()
    analyzer.estimate_pKa()
    analyzer.analyze_hydrogen_bonds()
    analyzer.decompose_dipole()
    
    # Or run everything at once:
    analyzer.analyze_all()
    
    # For molecules with defined groups:
    groups = {
        'methyl': [0, 1, 2, 3],      # CH3 group
        'carboxyl': [4, 5, 6, 7]     # COOH group
    }
    analyzer.analyze_charge_transfer(groups=groups)
    
    # Create charge flow animation (needs normal mode from Gaussian):
    # normal_mode = np.array([...])  # From Gaussian output
    # analyzer.animate_charge_flow(normal_mode, espaloma_calc)
    """
    
    print(example_code)
    
    print("\n" + "="*70)
    print("To integrate with your workflow:")
    print("  1. Import this module: from charge_analysis import ChargeAnalyzer")
    print("  2. Pass atoms and charges from Espaloma")
    print("  3. Call analysis functions as needed")
    print("="*70 + "\n")
