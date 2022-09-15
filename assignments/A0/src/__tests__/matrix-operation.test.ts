import { apply as _, pipe } from 'fp-ts/lib/function'
import { describe, expect, it } from 'vitest'
import * as M from '../matrixe-operation'

describe('multiply', () => {
  it('2x3 · 3x2', () => {
    expect(pipe(
      [
        [1, 1, 1],
        [1, 1, 1],
      ],
      M.multiply,
      _([
        [2, 2],
        [2, 2],
        [2, 2],
      ]),
    )).toEqual(
      [
        [6, 6],
        [6, 6],
      ],
    )
  })
})

describe('transpose', () => {
  it('2x3', () => {
    expect(
      M.transpose(
        [
          [1, 2, 3],
          [4, 5, 6],
        ],
      ),
    ).toEqual(
      [
        [1, 4],
        [2, 5],
        [3, 6],
      ],
    )
  })
})
