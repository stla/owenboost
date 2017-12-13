{-# LANGUAGE ForeignFunctionInterface #-}
module OwenExport
  where
import           Foreign
import           Foreign.C
import Test

foreign export ccall pnormExport :: Ptr CDouble -> IO ()
pnormExport :: Ptr CDouble -> IO ()
pnormExport result = do
  poke result pnorm

foreign export ccall studentCexport :: Ptr CDouble -> IO ()
studentCexport :: Ptr CDouble -> IO ()
studentCexport result = do
  (>>=) (studentC 1 1 [1,2]) (pokeArray result)
