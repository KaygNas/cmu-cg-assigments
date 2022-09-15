import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import * as N from '../../number-operation'
import type { Vector } from '../../interface'
import * as V from '../../vector-operation'

describe('exercise-1', () => {
  const u: Vector = [6, 2]
  const v: Vector = [8, 9]
  const a = 9
  const b = 3
  it('a', () => {
    expect(pipe(
      u, V.add, _(v),
    )).toEqual(
      [14, 11],
    )
  })
  it('b', () => {
    expect(pipe(
      b, V.scale, _(u),
    )).toEqual(
      [18, 6],
    )
  })
  it('c', () => {
    expect(pipe(
      V.scale(a)(u), V.subtract, _(V.scale(b)(v)),
    )).toEqual(
      [30, -9],
    )
  })
})

describe('exercise-2', () => {
  const u: Vector = [2, 2, 8]
  const v: Vector = [8, 8, 7]
  it('a', () => {
    expect(pipe(
      u, V.subtract, _(v),
    )).toEqual(
      [-6, -6, 1],
    )
  })
  it('b', () => {
    expect(pipe(
      u, V.add, _(V.scale(4)(v)),
    )).toEqual(
      [34, 34, 36],
    )
  })
})

describe('exercise-3', () => {
  const p = (x: number) => 2 * x ** 2 + 2 * x + 8
  const q = (x: number) => 8 * x ** 2 + 8 * x + 7
  const x = 4

  it('a', () => {
    expect(pipe(
      p(x), N.subtract, _(q(x)),
    )).toBe(
      -119,
    )
  })

  it('b', () => {
    expect(pipe(
      p(x), N.add, _(N.scale(4)(q(x))),
    )).toBe(
      716,
    )
  })
})
