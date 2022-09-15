import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import type { Vector } from '../../interface'
import * as V from '../../vector-operation'
import * as N from '../../number-operation'

describe('exercies-18', () => {
  const x: Vector = [9, 9]
  const y: Vector = [2, 3]
  const z: Vector = [3, 9]
  const a = 5
  const b = 7

  it('a-b', () => {
    const ra = pipe(
      pipe(a, V.scale, _(x)),
      V.add,
      _(pipe(b, V.scale, _(y))),
      V.innerProduct,
      _(z),
    )
    const rb = pipe(
      a, N.scale, _(pipe(x, V.innerProduct, _(z))),
      N.add,
      _(pipe(b, N.scale, _(pipe(y, V.innerProduct, _(z))))),
    )
    expect(ra).toBe(rb)
  })

  it('c-d', () => {
    const rc = V.norm(pipe(a, V.scale, _(x))) ** 2
    const rd = a ** 2 * V.norm(x) ** 2
    expect(rc).toBeCloseTo(rd)
  })
})
