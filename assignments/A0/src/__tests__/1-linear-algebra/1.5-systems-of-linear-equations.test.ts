import { describe, expect, it } from 'vitest'
import type { Vector } from '../../interface'

describe('exercise-16', () => {
  const f1 = (x: number, y: number) => 5 * x + 6 * y
  const f2 = (x: number, y: number) => -6 * x + 5 * y
  it('', () => {
    const x = 22 / 61
    const y = 63 / 61
    expect(
      f1(x, y),
    ).toBeCloseTo(
      8,
    )
    expect(
      f2(x, y),
    ).toBeCloseTo(
      3,
    )
  })
})

describe('exercise-17', () => {
  const f1 = ([x, y, z]: Vector) => x + y + z
  const f2 = ([x, y, z]: Vector) => x - y + z
  const f3 = ([x, y, z]: Vector) => x ** 2 + y ** 2 + z ** 2
  const z = ((29 ** (1 / 2)) + 1) / 8
  const y = -1 / 4
  const x = 1 / 4 - z
  const v: Vector = [x, y, z]

  it('', () => {
    expect(
      f1(v),
    ).toBe(
      0,
    )

    expect(
      f2(v),
    ).toBe(
      1 / 2,
    )

    expect(
      f3(v),
    ).toBeCloseTo(
      1,
    )
  })
})
