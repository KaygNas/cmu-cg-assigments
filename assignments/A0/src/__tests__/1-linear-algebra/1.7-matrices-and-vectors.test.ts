import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import { round } from 'lodash'
import type { Matrix, Vector } from '../../interface'
import * as M from '../../matrixe-operation'

describe('exercise-21', () => {
  const f = ([u1, u2]: Vector) => [9 * u1 + 2 * u2, 7 * u1 + 4 * u2]
  const m: Matrix = [
    [9, 2],
    [7, 4],
  ]
  it('a', () => {
    const u: Vector = [2, 1]
    expect(
      M.multiply(m)(M.v2m(u)),
    ).toEqual(
      M.v2m(f(u)),
    )
  })

  it('b', () => {
    const x: Vector = [6, 2]
    expect(
      M.multiply(m)(M.v2m(x)),
    ).toEqual(
      M.v2m([54 + 4, 42 + 8]),
    )
  })
})

describe('exercise-22', () => {
  type R2 = [number, number]
  type R3 = [number, number, number]
  const f = ([x1, x2]: R2): R3 => [3 * x1, 4 * x2, x1 + x2]
  const g = ([x1, x2, x3]: R3): R3 => [2 * x2, 2 * x3, 2 * x1]
  const mA: Matrix = [
    [3, 0],
    [0, 4],
    [1, 1],
  ]
  const mB: Matrix = [
    [0, 2, 0],
    [0, 0, 2],
    [2, 0, 0],
  ]
  it('a-b-c', () => {
    const mBA = M.multiply(mB)(mA)
    const x: R2 = [2, 4]
    expect(pipe(
      mBA, M.multiply, _(M.v2m(x)),
    )).toEqual(pipe(
      x, f, g, M.v2m,
    ))
  })
})

describe('exercise-23', () => {
  const mI: Matrix = [
    [1, 0],
    [0, 1],
  ]
  it('a', () => {
    const mA: Matrix = [
      [2, 3],
      [4, 5],
    ]
    expect(pipe(
      mA, M.multiply, _(mI),
    )).toEqual(
      mA,
    )
  })

  const f = ([x, y]: Vector) => [3 * x + 3 * y, -6 * x + 9 * y]
  const a = 3
  const b = 3
  const c = -6
  const d = 9
  const mA: Matrix = [
    [a, b],
    [c, d],
  ]
  it('b', () => {
    const u: Vector = [3, 6]
    expect(pipe(
      mA, M.multiply, _(M.v2m(u)),
    )).toEqual(
      M.v2m(f(u)),
    )
  })

  it('c', () => {
    const bc减ad = b * c - a * d
    const mAI: Matrix = [
      [-d / bc减ad, b / bc减ad],
      [c / bc减ad, -a / bc减ad],
    ]
    expect(pipe(
      mA, M.multiply, _(mAI),
      M.traverseMatrix((m, n, v) => Math.ceil(v) || 0),
    )).toEqual(
      mI,
    )
  })
})

describe('exercise-24', () => {
  const mA: Matrix = [
    [3, 8],
    [-8, 3],
  ]
  const mAT: Matrix = [
    [3, -8],
    [8, 3],
  ]
  it('a', () => {
    expect(pipe(
      mA, M.transpose,
    )).toEqual(
      mAT,
    )
  })

  it('b', () => {
    /**
     * 两个相互垂直且长度为 a 的向量构成的矩阵，
     * 它们相乘的结果为 a * IdentityMatrix。
     * (IdentityMatrix 为单位元矩阵，即对任意矩阵 M, M * IdentityMatrix = M)
     */
    expect(pipe(
      mA, M.multiply, _(mAT),
    )).toEqual([
      [73, 0],
      [0, 73],
    ])
  })
})

describe('exercise-25', () => {
  it('a', () => {
    const baseUV = 97 ** (1 / 2)
    const mUV: Matrix = [
      [9 / baseUV, -4 / baseUV],
      [4 / baseUV, 9 / baseUV],
    ]

    const baseXY = 53 ** (1 / 2)
    const mXY: Matrix = [
      [7 / baseXY, -2 / baseXY],
      [2 / baseXY, 7 / baseXY],
    ]
    const E = pipe(
      mXY, M.multiply, _(M.transpose(mUV)),
    )

    const roundMatrix = M.traverseMatrix((m, n, v) => round(v, 12))
    expect(pipe(
      E, M.multiply, _(mUV),
      M.scale(10),
      M.traverseMatrix((m, n, v) => round(v, 12)),
      roundMatrix,
    )).toEqual(pipe(
      mXY,
      M.scale(10),
      roundMatrix,
    ),
    )
  })
})
