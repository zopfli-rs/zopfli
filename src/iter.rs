use std::iter::Peekable;

pub struct IsFinal<I>
where
    I: Iterator,
{
    inner: Peekable<I>,
}

impl<I> Iterator for IsFinal<I>
where
    I: Iterator,
{
    type Item = (I::Item, bool);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner
            .next()
            .map(|val| (val, self.inner.peek().is_none()))
    }
}

pub trait IsFinalIterator: Iterator + Sized {
    fn is_final(self) -> IsFinal<Self> {
        IsFinal {
            inner: self.peekable(),
        }
    }
}

impl<T> IsFinalIterator for T where T: Iterator + Sized {}
