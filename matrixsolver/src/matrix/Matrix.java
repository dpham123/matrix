package matrix;

import java.util.*;

public class Matrix{
	private double[][] matrix;

	public Matrix(double[][] matrix) throws InvalidMatrixException {
		for (int i = 0; i < matrix.length; i++) {
			if (matrix[0].length != matrix[i].length) {
				throw new InvalidMatrixException("Invalid matrix. Saw " + matrix[0].length + " columns on row 0 and "
						+ matrix[i].length + " columns on row " + i);
			}
		}
		this.matrix = matrix; // Arrays whose slots have nothing are initialized to 0.0
	}
	/**
	 * 
	 * @return matrix with no answer column
	 * @throws InvalidMatrixException
	 */
	private Matrix getMatrixNoAns() throws InvalidMatrixException {
		double[][] matrix = new double[this.matrix.length][this.matrix[0].length - 1];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				matrix[i][j] = this.matrix[i][j];
			}
		}
		return new Matrix(matrix);
	}

	public String toString() {
		String s = "";
		for (int row = 0; row < matrix.length; row++) {
			for (int column = 0; column < matrix[0].length; column++) {
				s += (matrix[row][column] + " ");
			}
			s += "\n";
		}
		return s;
	}

	private double getDeterminant() throws DeterminantException {
		// Checks if matrix is an n x n matrix
		for (int i = 0; i < this.matrix.length; i++) {
			if (this.matrix.length != this.matrix[i].length) {
				throw new DeterminantException("Matrix does not have the same # of rows as columns");
			}
		}
		return recurseDeterminant(this.matrix);
	}

	private double recurseDeterminant(double[][] matrix) {
		double[] coefficients = new double[matrix.length];
		double[][][] subMatrices = new double[matrix.length][matrix.length][matrix.length];

		// Checks if matrix is 2x2
		if (matrix.length == 2) {
			return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
		} else {
			for (int i = 0; i < matrix.length; i++) {
				coefficients[i] = matrix[0][i];
				subMatrices[i] = subMatrix(i, matrix);
			}

			double determinant = 0;
			int counter = 1;
			for (int i = 0; i < coefficients.length; i++) {
				// Checks for place of sign
				if (counter % 2 == 1) {
					determinant += coefficients[i] * recurseDeterminant(subMatrices[i]);
				} else {
					determinant += -coefficients[i] * recurseDeterminant(subMatrices[i]);
				}
				counter++;
			}
			return determinant;
		}
	}

	private double[][] cramerSubMatrix(int columnToReplace, double[][] matrix, double[] rightHandNumbers) {
		// Populates submatrix array with numbers from matrix
		double[][] subMatrix = new double[matrix.length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				subMatrix[i][j] = matrix[i][j];
			}
		}

		// Replaces indicated column with right hand numbers
		for (int i = 0; i < subMatrix.length; i++) {
			subMatrix[i][columnToReplace] = rightHandNumbers[i];
		}
		return subMatrix;
	}

	private static double[][] subMatrix(int columnToSkip, double[][] matrix) {
		double[][] subMatrix = new double[matrix.length - 1][matrix.length - 1];
		double[] valuesOfSubMatrix = new double[(matrix.length - 1) * (matrix.length - 1)];

		// Looks through matrix to find numbers to add to submatrix
		int counter = 0;
		for (int i = 1; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				if (j != columnToSkip) {
					valuesOfSubMatrix[counter] = matrix[i][j];
					counter++;
				}
			}
		}

		// Populates subMatrix with stored matrix numbers
		counter = 0;
		for (int i = 0; i < subMatrix.length; i++) {
			for (int j = 0; j < subMatrix.length; j++) {
				subMatrix[i][j] = valuesOfSubMatrix[counter];
				counter++;
			}
		}
		return subMatrix;
	}

	// Solves matrix using Cramer's Rule
	private ArrayList<Double> solveMatrixCramer() throws UnsolvableMatrixException {
		if (this.matrix.length != this.matrix[0].length - 1) {
			throw new UnsolvableMatrixException("Expected n x (n + 1) matrix. Saw " + this.matrix.length + " x "
					+ this.matrix[0].length + " matrix");
		}

		// Constructs matrix without right hand numbers to calculate determinants
		double[][] subMatrix = new double[this.matrix.length][this.matrix[0].length - 1];
		for (int i = 0; i < this.matrix.length; i++) {
			for (int j = 0; j < this.matrix[0].length - 1; j++) {
				subMatrix[i][j] = this.matrix[i][j];
			}
		}

		// Checks whether matrix is solvable
		double determinant = recurseDeterminant(subMatrix);
		if (determinant == 0) {
			throw new UnsolvableMatrixException("Matrix is not invertible and not solvable");
		}

		// Populates array with right hand numbers
		double[] rightHandNumbers = new double[this.matrix.length];
		for (int i = 0; i < this.matrix.length; i++) {
			rightHandNumbers[i] = this.matrix[i][this.matrix.length];
		}

		// Creates "Cramer" submatrices
		double[][][] subMatrices = new double[subMatrix[0].length][subMatrix[0].length][subMatrix[0].length];
		for (int j = 0; j < subMatrix[0].length; j++) {
			subMatrices[j] = cramerSubMatrix(j, subMatrix, rightHandNumbers);
		}

		// Computes determinants for individual "Cramer" submatrices and returns answer
		ArrayList<Double> solutions = new ArrayList<>();
		for (double[][] d : subMatrices) {
			solutions.add(recurseDeterminant(d) / determinant);
		}
		return solutions;
	}

	private ArrayList<Double> solveMatrixGaussian() throws UnsolvableMatrixException {
		if (this.matrix.length != this.matrix[0].length - 1) {
			throw new UnsolvableMatrixException("Expected n x (n + 1) matrix. Saw " + this.matrix.length + " x "
					+ this.matrix[0].length + " matrix");
		}

		for (int i = 0; i < matrix.length; i++) {
			clearColumn(i);
		}

		// Populates arraylist with right hand numbers
		ArrayList<Double> rightHandNumbers = new ArrayList<>();
		for (int i = 0; i < this.matrix.length; i++) {
			rightHandNumbers.add(this.matrix[i][this.matrix.length]);
		}
		return rightHandNumbers;
	}

	private void clearColumn(int column) throws UnsolvableMatrixException {
		// Switches rows if the number cannot be made into a leading 1
		if (this.matrix[column][column] == 0) {
			for (int i = column + 1; i < this.matrix.length; i++) {
				if (this.matrix[i][column] != 0) {
					rowSwitch(i, column);
				}
			}
		}

		// Gets leading 1 in the column
		double numberToDivideBy = this.matrix[column][column];
		if (numberToDivideBy == 0) {
			throw new UnsolvableMatrixException("Matrix has no solution");
		}
		if (numberToDivideBy != 1) {
			for (int j = 0; j < this.matrix[0].length; j++) {
				this.matrix[column][j] /= numberToDivideBy;
			}
		}

		// Clears out the rest of the column with 0's
		double coefficient = -1;
		boolean coefficientDetermined = false;
		for (int i = 0; i < this.matrix.length; i++) {
			if (this.matrix[i][column] == 0 || i == column) {
				continue;
			}
			
			for (int j = 0; j < this.matrix[0].length; j++) {
				if (!coefficientDetermined) {
					coefficient = this.matrix[i][column] * -this.matrix[column][column];
					coefficientDetermined = true;
				}
				this.matrix[i][j] += coefficient * this.matrix[column][j];
			}
			coefficientDetermined = false;
		}
	}

	private void rowSwitch(int firstRow, int secondRow) {
		double[] tempRow = this.matrix[firstRow];
		this.matrix[firstRow] = this.matrix[secondRow];
		this.matrix[secondRow] = tempRow;
	}

	private static void sop(Object x) {
		System.out.println(x);
	}

	private static int getRandomNumberInRange(int min, int max) {

		if (min >= max) {
			throw new IllegalArgumentException("Max must be greater than min");
		}

		Random r = new Random();
		return r.nextInt((max - min) + 1) + min;
	}

	public static void main(String[] args) {

		double[][] testMatrix1 = { { 37, 15, 135, 25, 20, 4, 6, 3, 6, 3 }, { 100, 44, 73, 84, 42, 4, 7, 4, 5, 21 },
				{ 877, 43, 47, 78, 34, 4, 6, 23, 56, 4 }, { 87, 74, 37, 98, 44, 3, 7, 5, 34, 6 },
				{ 875, 94, 75, 98, 54, 3, 5, 7, 3, 4 }, { 4, 96, 55, 95, 524, 4, 2, 6, 8, 2 },
				{ 64, 1, 75, 58, 84, 4, 7, 2, 8, 6 }, { 84, 91, 75, 98, 4, 3, 5, 7, 6, 3 },
				{ 86, 34, 75, 98, 3, 3, 3, 2, 37, 9 } };

		double[][] v1 = { { 0, 1, 8, 6, 4 }, { 0, 6, 7, 4, 8 }, { 0, 4, 3, 2, 1 }, { 4, 3, 8, 9, 5 } };
		double[][] v3 = { { 3, 0, 1, 5 }, { 2, 7, 6, 9 }, { 3, 6, 1, 3 } };
		double[][] v4 = { { 0, -1, 1 }, { 0, -1, 1 }, { 1, -2, 0 } };
		
		try {
			new Matrix(v1);
		} catch (InvalidMatrixException e3) {
			e3.printStackTrace();
		}
		for (int trials = 0; trials < 100000; trials++) {
			int a = getRandomNumberInRange(0, 9);
			int b = getRandomNumberInRange(0, 9);
			int c = getRandomNumberInRange(0, 9);
			int d = getRandomNumberInRange(0, 9);
			int e = getRandomNumberInRange(0, 9);
			int f = getRandomNumberInRange(0, 9);
			int g = getRandomNumberInRange(0, 9);
			int h = getRandomNumberInRange(0, 9);
			int i = getRandomNumberInRange(0, 9);
			int j = getRandomNumberInRange(0, 9);
			int k = getRandomNumberInRange(0, 9);
			int l = getRandomNumberInRange(0, 9);
			int m = getRandomNumberInRange(0, 9);
			int n = getRandomNumberInRange(0, 9);
			int o = getRandomNumberInRange(0, 9);
			int p = getRandomNumberInRange(0, 9);

			double[][] v2 = { { a, b, c, d, 5 }, { e, f, g, h, 9 }, { i, j, k, l, 3 }, { m, n, o, p, 7 }};
			double[][] v5 = { { a, b, c, d, 5 }, { e, f, g, h, 9 }, { i, j, k, l, 3 }, { m, n, o, p, 7 }};
			try {
				// Initializes two identical matrices
				Matrix matrix2 = new Matrix(v2);
				Matrix original = new Matrix(v5);

				// Tries to solve one of the matrices
				try {
					matrix2.solveMatrixGaussian();
				} catch (UnsolvableMatrixException e1) {

					// Gets the determinant of the failed matrix
					try {
						double determinant = original.getMatrixNoAns().getDeterminant();
						if (determinant != 0) {
							sop("Original: ");
							sop(original);
							sop("Final: ");
							sop(matrix2);
							sop("Determinant: ");
							sop(determinant);
							System.out.println();
						}
					} catch (DeterminantException e2) {
						e2.printStackTrace();
					}
				}
			} catch (InvalidMatrixException e1) {
				e1.printStackTrace();
			}
		}
		sop("No bugs :)");
	}

}
