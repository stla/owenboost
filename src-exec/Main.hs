module Main
  where
import Test

main :: IO()
main = do
  x <- studentC 1 2 [1,2,3]
  print x
