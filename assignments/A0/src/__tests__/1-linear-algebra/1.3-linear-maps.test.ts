import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import type { Vector } from '../../interface'
import * as V from '../../vector-operation'

describe('exercise-9', () => {
  /**
   * not linear
   */
  const f = (x: number) => 2 * x + 2
  /**
   * linear
   */
  const g = (x: number) => 7 * x
  /**
   * not linear
   */
  const h = (x: number) => -1 * x ** 2
  const x = 7
  const y = 9

  it('a', () => {
    expect(
      f(x + y) - (f(x) + f(y)),
    ).not.toBe(
      0,
    )
  })

  it('b', () => {
    expect(
      f(6 * x) - 6 * f(x),
    ).not.toBe(
      0,
    )
  })

  it('c', () => {
    expect(g(x + y) - (g(x) + g(y))).toBe(0)
  })

  it('d', () => {
    expect(g(5 * x) - 5 * g(x)).toBe(0)
  })

  it('e', () => {
    expect(h(x + y) - h(x) + h(y)).not.toBe(0)
  })

  it('f', () => {
    expect(h(2 * x) - 2 * h(x)).not.toBe(0)
  })
})

describe('exercies-10', () => {
  const f = ([u1, u2]: Vector) => u1 + u2 + 8
  const g = ([u1, u2]: Vector) => u1 ** 2 + u2 ** 2
  const u: Vector = [7, 7]
  const v: Vector = [3, 8]
  const w1 = 0.3
  const w2 = 0.7
  it('a', () => {
    expect(f(
      pipe(
        w1, V.scale, _(u),
        V.add,
        _(pipe(
          w2, V.scale, _(v),
        )),
      ))
      - (w1 * f(u) + w2 * f(v)),
    ).toBe(
      0,
    )
  })

  it('b', () => {
    expect(g(
      pipe(
        w1, V.scale, _(u),
        V.add,
        _(pipe(
          w2, V.scale, _(v),
        )),
      )) - (w1 * g(u) + w2 * g(v)),
    ).not.toBe(
      0,
    )
  })
})
