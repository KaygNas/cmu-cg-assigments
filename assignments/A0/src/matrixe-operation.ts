import { apply, pipe } from 'fp-ts/lib/function'
import type { Matrix, Vector } from './interface'
import * as V from './vector-operation'

const rowLength = (m: Matrix): number => {
	return m.length
}
const colLength = (m: Matrix): number => {
	return m[0].length
}

const nthRowOf = (n: number, m: Matrix): Vector => {
	return m[n]
}
const nthColOf = (n: number, m: Matrix): Vector => {
	return m.map((row) => row[n])
}

const makeEmptyMatrix = (m: number, n: number): Matrix => {
	return Array.from({ length: m }).map(() => Array.from({ length: n }))
}

type Iteratee = (m: number, n: number, v: number) => number
export const traverseMatrix = (iteratee: Iteratee) => (matrix: Matrix) => {
	return matrix.map((row, m) => {
		return row.map((v, n) => {
			return iteratee(m, n, v)
		})
	})
}

export const multiply = (m1: Matrix) => (m2: Matrix) => {
	return traverseMatrix((m, n) => {
		return pipe(nthRowOf(m, m1), V.innerProduct, apply(nthColOf(n, m2)))
	})(makeEmptyMatrix(rowLength(m1), colLength(m2)))
}

export const transpose = (matrix: Matrix): Matrix => {
	return traverseMatrix((m, n) => matrix[n][m])(
		makeEmptyMatrix(colLength(matrix), rowLength(matrix)),
	)
}

export const v2m = (v: Vector): Matrix => {
	return transpose([v])
}

export const scale = (scalar: number) => (m: Matrix) => {
	return traverseMatrix((m, n, v) => scalar * v)(m)
}
