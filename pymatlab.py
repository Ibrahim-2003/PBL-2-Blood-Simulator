import pyperclip as pc
texts = ''
for i in range(1,61):
    texts += f'eqn{i} '

pc.copy(texts)
print(texts)