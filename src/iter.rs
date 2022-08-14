use std::iter::Peekable;

pub struct FlagLastIterator<I>
where
    I: Iterator,
{
    inner: Peekable<I>,
}

impl<I> Iterator for FlagLastIterator<I>
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

pub trait ToFlagLastIterator: Iterator + Sized {
    fn flag_last(self) -> FlagLastIterator<Self> {
        FlagLastIterator {
            inner: self.peekable(),
        }
    }
}

impl<T> ToFlagLastIterator for T where T: Iterator + Sized {}
