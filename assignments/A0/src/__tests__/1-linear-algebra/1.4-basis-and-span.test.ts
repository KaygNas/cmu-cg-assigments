import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import type { Vector } from '../../interface'
import * as V from '../../vector-operation'

describe('exercise-12', () => {
  const e1: Vector = [1 / (2 ** (1 / 2)), 1 / (2 ** (1 / 2))]
  const e2: Vector = [-1 / (2 ** (1 / 2)), 1 / (2 ** (1 / 2))]
  const u: Vector = [3, 6]

  const a = pipe(
    e1, V.innerProduct, _(u),
  )
  it('a', () => {
    expect(a).toBe(
      9 / (2 ** (1 / 2)),
    )
  })

  const b = pipe(
    e2, V.innerProduct, _(u),
  )
  it('b', () => {
    expect(b).toBe(
      3 / (2 ** (1 / 2)),
    )
  })

  it('c', () => {
    const result = pipe(
      u, V.subtract, _(pipe(
        pipe(a, V.scale, _(e1)),
        V.add,
        _(pipe(b, V.scale, _(e2))),
      )),
    )
    // 存在浮点计算精度问题，只能取近似值
    expect(
      result[0],
    ).toBeCloseTo(
      0,
    )
    expect(
      result[1],
    ).toBeCloseTo(
      0,
    )
  })
})

describe('exercies-13', () => {
  it('', () => {
    const e1: Vector = [1, 0]
    const e2: Vector = [1, 1]
    const u: Vector = [3, 3]
    const a = pipe(
      e1, V.innerProduct, _(u),
    )
    const b = pipe(
      e2, V.innerProduct, _(u),

    )
    const result = pipe(
      pipe(a, V.scale, _(e1)),
      V.add,
      _(pipe(b, V.scale, _(e2))),
    )

    expect(
      result[0],
    ).toBeCloseTo(
      9,
    )
    expect(
      result[1],
    ).toBeCloseTo(
      6,
    )
  })
})

describe('exercise-15', () => {
  const e1: Vector = [4, 9]
  const e2: Vector = [8, 18]
  const w: Vector = [2, 7]
  const det = ([u1, u2]: Vector, [v1, v2]: Vector) => u1 * v2 - u2 * v1
  it('a', () => {
    expect(
      det(e1, e2),
    ).toBe(
      // 4 * 18 - 9 * 8,
      0,
    )
  })

  const a = pipe(
    e1, V.innerProduct, _(w),
  )
  it('b', () => {
    expect(
      a,
    ).toBe(
      71,
    )
  })

  const b = pipe(
    e2, V.innerProduct, _(w),
  )
  it('c', () => {
    expect(
      b,
    ).toBe(
      142,
    )
  })

  it('d', () => {
    expect(pipe(
      a, V.scale, _(e1),
      V.add,
      _(pipe(b, V.scale, _(e2))),
    )).toEqual(
      [1420, 3195],
    )
  })
})
