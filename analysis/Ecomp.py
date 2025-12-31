import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from statsmodels.stats.multitest import multipletests
import pingouin as pg
import sys
import argparse
import logging
from itertools import combinations
from pathlib import Path
from typing import Optional, Dict
import json


class EcompAnalyzer:
    """Analyzer for E_comp data with statistical analysis and visualization."""

    def __init__(self, 
                 file_path: str, 
                 alpha: float = 0.05, 
                 test_normality: bool = True, 
                 export_path: Optional[str] = None,
                 outlier_threshold: float = 3.0,
                 save_plots: bool = True,
                 detect_outliers: bool = True):
        """
        Initialize EcompAnalyzer.

        Args:
            file_path: Path to input CSV file
            alpha: Significance level
            test_normality: Whether to perform normality tests
            export_path: Path for exporting results
            outlier_threshold: Z-score threshold for outlier detection
            save_plots: Whether to save plots
            detect_outliers: Whether to detect outliers
        """
        self.file_path = Path(file_path)
        self.export_path = Path(export_path) if export_path else self.file_path.parent
        self.alpha = alpha
        self.test_normality = test_normality
        self.outlier_threshold = outlier_threshold
        self.save_plots = save_plots
        self.detect_outliers = detect_outliers
        self.data = None
        self.results: Dict = {}
        
        # Set plotting style
        # plt.style.use('seaborn')  # optional; depends on matplotlib version
        # Configure logging
        self._setup_logging()

    def _setup_logging(self) -> None:
        """Configure logging settings."""
        log_file = self.export_path / 'analysis.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )

    def load_and_validate_data(self) -> pd.DataFrame:
        """Load and validate input data."""
        try:
            self.data = pd.read_csv(self.file_path)
            self._validate_columns()
            self._validate_numeric_data()
            self._calculate_ecomp()
            if self.detect_outliers:
                self._detect_outliers()
            return self.data
        except Exception as e:
            logging.error(f"Error loading data: {str(e)}")
            raise

    def _validate_columns(self) -> None:
        """Validate required columns exist."""
        required_columns = ['E1', 'E2', 'Group']
        if not all(col in self.data.columns for col in required_columns):
            raise ValueError(f"Missing required columns: {required_columns}")

    def _validate_numeric_data(self) -> None:
        """Validate numeric columns contain valid data."""
        numeric_cols = ['E1', 'E2']
        for col in numeric_cols:
            if not pd.to_numeric(self.data[col], errors='coerce').notnull().all():
                raise ValueError(f"Column {col} contains non-numeric values")
            if (self.data[col] <= 0).any():
                raise ValueError(f"Column {col} contains non-positive values")

    def _calculate_ecomp(self) -> None:
        """Calculate E_comp and standardized version."""
        self.data['E_comp'] = (2 * self.data['E1'] * self.data['E2']) / (self.data['E1'] + self.data['E2'])
        self.data['E_comp_z'] = stats.zscore(self.data['E_comp'])

    def _detect_outliers(self) -> None:
        """Detect outliers using z-score method."""
        z_scores = np.abs(stats.zscore(self.data['E_comp']))
        self.data['is_outlier'] = z_scores > self.outlier_threshold
        n_outliers = sum(self.data['is_outlier'])
        logging.info(f"Data loaded: {len(self.data)} samples, {n_outliers} outliers")

    def plot_distributions(self) -> None:
        """Generate distribution plots."""
        try:
            fig = plt.figure(figsize=(15, 12))
            gs = fig.add_gridspec(3, 3)

            self._plot_boxplot(fig.add_subplot(gs[0, 0]))
            self._plot_violin(fig.add_subplot(gs[0, 1]))
            self._plot_kde(fig.add_subplot(gs[0, 2]))
            self._plot_qq(fig.add_subplot(gs[1, :]))
            self._plot_residuals(fig.add_subplot(gs[2, :]))

            plt.tight_layout()

            if self.save_plots:
                plt.savefig(self.export_path / 'distributions.png', dpi=300, bbox_inches='tight')
            plt.show()

        except Exception as e:
            logging.error(f"Error in plotting: {str(e)}")
            raise

    def _plot_boxplot(self, ax: plt.Axes) -> None:
        sns.boxplot(data=self.data, x='Group', y='E_comp', ax=ax)
        ax.set_title('Distribution by Group')
        ax.tick_params(axis='x', rotation=45)

    def _plot_violin(self, ax: plt.Axes) -> None:
        sns.violinplot(data=self.data, x='Group', y='E_comp', ax=ax)
        ax.set_title('Violin Plot')
        ax.tick_params(axis='x', rotation=45)

    def _plot_kde(self, ax: plt.Axes) -> None:
        sns.kdeplot(data=self.data, x='E_comp', hue='Group', fill=True, ax=ax)
        ax.set_title('Density Distribution')

    def _plot_qq(self, ax: plt.Axes) -> None:
        for group in self.data['Group'].unique():
            group_data = self.data[self.data['Group'] == group]['E_comp']
            stats.probplot(group_data, dist="norm", plot=ax)
        ax.set_title('Q-Q Plots by Group')

    def _plot_residuals(self, ax: plt.Axes) -> None:
        residuals = self.data['E_comp'] - self.data.groupby('Group')['E_comp'].transform('mean')
        sns.scatterplot(data=self.data, x='E_comp', y=residuals, hue='Group', ax=ax)
        ax.axhline(y=0, color='r', linestyle='--')
        ax.set_title('Residual Plot')
        if self.save_plots:
            plt.savefig(self.export_path / 'residuals.png', dpi=300, bbox_inches='tight')

    def perform_statistical_analysis(self) -> Dict:
        """Perform statistical tests and post hoc analyses."""
        try:
            groups = self.data['Group'].unique()
            group_data = [self.data[self.data['Group'] == group]['E_comp'] for group in groups]

            # Normality tests
            if self.test_normality:
                self.results['normality'] = {
                    group: stats.shapiro(data) for group, data in zip(groups, group_data)
                }
                is_normal = all(p > self.alpha for _, (_, p) in self.results['normality'].items())
            else:
                is_normal = True

            # Homogeneity of variance
            _, levene_p = stats.levene(*group_data)
            self.results['homogeneity'] = {'statistic': _, 'p_value': levene_p}

            # Main statistical test
            if is_normal and levene_p > self.alpha:
                f_stat, p_value = stats.f_oneway(*group_data)
                test_used = "ANOVA"
                effect_size = pg.anova(data=self.data, dv='E_comp', between='Group')['np2'][0]

                # Tukey's HSD post hoc test
                post_hoc = pg.pairwise_tukey(data=self.data, dv='E_comp', between='Group')
            else:
                h_stat, p_value = stats.kruskal(*group_data)
                test_used = "Kruskal-Wallis"
                effect_size = pg.kruskal(data=self.data, dv='E_comp', between='Group')['np2'][0]

                # Dunn's Test post hoc
                post_hoc = pg.pairwise_dunn(data=self.data, dv='E_comp', between='Group', p_adjust='bonferroni')

            # Store main test results
            self.results.update({
                'test_used': test_used,
                'main_test': {'statistic': f_stat if test_used == "ANOVA" else h_stat, 'p_value': p_value},
                'effect_size': effect_size,
                'pairwise_tests': post_hoc.to_dict('records')
            })

            # Save post hoc heatmap
            self.plot_post_hoc_matrix(post_hoc)
            return self.results
        except Exception as e:
            logging.error(f"Error in statistical analysis: {str(e)}")
            raise

    def plot_post_hoc_matrix(self, post_hoc: pd.DataFrame) -> None:
        """Generate and save a heatmap of post hoc corrected p-values."""
        try:
            matrix = post_hoc.pivot(index='A', columns='B', values='p-corr')
            plt.figure(figsize=(10, 8))
            sns.heatmap(matrix, annot=True, fmt=".3f", cmap='coolwarm', cbar_kws={'label': 'Corrected p-value'})
            plt.title("Pairwise Comparison Heatmap (Corrected p-values)")
            if self.save_plots:
                plt.savefig(self.export_path / 'post_hoc_matrix.png', dpi=300, bbox_inches='tight')
            plt.show()
        except Exception as e:
            logging.error(f"Error plotting post hoc heatmap: {str(e)}")
            raise

    def export_results(self) -> None:
        """Export results to JSON and CSV."""
        try:
            # Save detailed results as JSON
            with open(self.export_path / 'analysis_results.json', 'w') as f:
                json.dump(self.results, f, indent=4, default=str)

            # Save summary as CSV
            summary = pd.DataFrame({
                'Test': [self.results['test_used']],
                'Statistic': [self.results['main_test']['statistic']],
                'P-value': [self.results['main_test']['p_value']],
                'Effect Size': [self.results['effect_size']]
            })
            summary.to_csv(self.export_path / 'analysis_summary.csv', index=False)

            logging.info(f"Results exported to {self.export_path}")

        except Exception as e:
            logging.error(f"Error exporting results: {str(e)}")
            raise

    def run_analysis(self) -> Dict:
        """Run complete analysis pipeline."""
        try:
            self.load_and_validate_data()
            self.plot_distributions()
            results = self.perform_statistical_analysis()
            self.export_results()
            self._print_summary(results)
            return results
        except Exception as e:
            logging.error(f"Analysis failed: {str(e)}")
            raise

    def _print_summary(self, results: Dict) -> None:
        """Print analysis summary."""
        print("\nAnalysis Summary:")
        print(f"Test used: {results['test_used']}")
        print(f"Statistic: {results['main_test']['statistic']:.4f}")
        print(f"P-value: {results['main_test']['p_value']:.4f}")
        print(f"Effect Size: {results['effect_size']:.4f}")
        print("\nPairwise Comparisons:")
        for test in results['pairwise_tests']:
            print(f"{test['A']} vs {test['B']}:")
            print(f"  Corrected p-value: {test['p-corr']:.4f}")
            print(f"  Effect size: {test['hedges']:.4f}")


def main():
    parser = argparse.ArgumentParser(description="E_comp analysis from CSV with columns E1,E2,Group")
    parser.add_argument("csv", help="Input CSV path")
    parser.add_argument("--out", default=None, help="Output directory for results (defaults to CSV directory)")
    args = parser.parse_args()

    analyzer = EcompAnalyzer(args.csv, export_path=args.out)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()