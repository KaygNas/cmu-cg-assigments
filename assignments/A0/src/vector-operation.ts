import type { Vector } from './interface'

export const scale = (a: number) => (u: Vector) => {
	return u.map((ui) => a * ui)
}

export const add = (u: Vector) => (v: Vector) => {
	return u.map((ui, i) => ui + v[i])
}

export const subtract = (u: Vector) => (v: Vector) => {
	return add(u)(scale(-1)(v))
}

export const innerProduct = (u: Vector) => (v: Vector) => {
	return u.reduce((result, ui, i) => result + ui * v[i], 0)
}

export const norm = (v: Vector) => {
	return innerProduct(v)(v) ** (1 / 2)
}
