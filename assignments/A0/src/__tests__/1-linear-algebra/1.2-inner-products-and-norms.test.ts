import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import type { Vector } from '../../interface'
import * as V from '../../vector-operation'
import * as N from '../../number-operation'

describe('exercise-4', () => {
  it('', () => {
    const u: Vector = [2, 5, 3]
    const v: Vector = [9, 4, 5]
    expect(pipe(
      u, V.innerProduct, _(v),
    )).toBe(
      18 + 20 + 15,
    )
  })
})

describe('exercise-5', () => {
  it('', () => {
    const v: Vector = [4, 6, 2]
    expect(pipe(
      V.norm(v),
    )).toBe(
      56 ** (1 / 2),
    )
  })
})

describe('exercise-6', () => {
  const ao = ([u1, u2]: Vector) => ([v1, v2]: Vector) => 6 * u1 * v1 + u1 * v2 + u2 * v1 + 5 * u2 * v2
  it('a: <x, x>', () => {
    /**
     * <x, x> >= 0
     */
    const x: Vector = [1, 0]
    expect(pipe(
      x, ao, _(x),
    )).toBe(
      6,
    )
  })

  it('b: <y, y>', () => {
    /**
     * <y, y> >= 0
     */
    const y: Vector = [0, 1]
    expect(pipe(
      y, ao, _(y),
    )).toBe(
      5,
    )
  })

  it('c: <u,v> - <v, u>', () => {
    /**
     * <u,v> = <v,u>
     */
    const u: Vector = [3, 9]
    const v: Vector = [5, 2]
    expect(pipe(
      pipe(u, ao, _(v)),
      N.subtract,
      _(pipe(v, ao, _(u))),
    )).toBe(
      0,
    )
  })

  it('d: <5u+v, w> - (5<u, w> + <v, w>)', () => {
    /**
     * <5u+v, w> - (5<u, w> + <v, w>)
     * = <5u,w> + <v, w> - (5<u, w> + <v, w>)
     * = 5<u, w> + <v, w> - (5<u, w> + <v, w>)
     * = 0
     */
    const u: Vector = [5, 4]
    const v: Vector = [5, 7]
    const w: Vector = [8, 4]

    expect(pipe(
      pipe(5, V.scale, _(u), V.add, _(v)),
      ao,
      _(w),
      N.subtract,
      _(pipe(
        5,
        N.scale,
        _(pipe(u, ao, _(w))),
        N.add,
        _(pipe(v, ao, _(w))),
      )),
    )).toBe(
      0,
    )
  })
})

describe('exercise-7', () => {
  it('a', () => {
    // TODO
  })
})

describe('exercise-8', () => {
  it('', () => {
    /**
    * f(x) = 6e^8x
    * 复合函数求导:
    * f(x) = g(h(x))
    * f'(x) = g'(h(x)) * h'(x)
    *
    * 设 g(x) = e^x, h(x) = 16x, 可得:
    * f(x) 的积分方程为: 36 * e^15x
    */
  })
})
